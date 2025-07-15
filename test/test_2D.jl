@testset "Params" begin
    params = TDGL_Polycrystal.get_params(tstep,GL)
    @test isa(params,Parameters)
    @test typeof(params.k) == Float64
    @test params.k == tstep
end
const params = TDGL_Polycrystal.get_params(tstep,GL)

@testset "Mesh" begin
    mesh = TDGL_Polycrystal.get_mesh(pixels,xmin,ymin,false)
    @test isa(mesh,RectMesh)
    @test isa(mesh.h,Tuple{Float64,Float64})
    @test mesh.extent == (xmin,ymin).*pixels
    @test mesh.periodic == (true,false)
end
const mesh = TDGL_Polycrystal.get_mesh(pixels,xmin,ymin,yperiodic)

const n0 = elemcount(mesh, Val(0))
const n1 = elemcount(mesh, Val(1))
const n2 = elemcount(mesh, Val(2))

@testset "State" begin
    state = TDGL_Polycrystal.get_state(n0,n1,mesh,backend)
    @test isa(state,State)
    @test isa(state.φ,RectPrimalForm0Data)
    @test sum(state.a.data) == 0
end
const state = TDGL_Polycrystal.get_state(n0,n1,mesh,backend)

@testset "Material" begin
    np_mesh = TDGL_Polycrystal.get_mesh(pixels,xmin,ymin,false)
    material,init_α,start_m⁻¹,start_σ = TDGL_Polycrystal.get_material(1.0,1.0,1.0,init_σ,n0,n1,pixels,grain_size,crystalangle,3,norm_res,norm_mass,alphaN,1.0,np_mesh,backend)

    @test isa(material,Material)
    @test isa(material.σ,RectPrimalForm1Data)
    @test isa(material.β,RectPrimalForm0Data)
    @test asarray(Adapt.adapt_structure(CPU(),material.α))[1,1] == -1.0

    @test length(material.σ.data) == elemcount(np_mesh, Val(1))
    @test sum(init_α)/n0 < 1
end
const material,α,β,m⁻¹,σ = TDGL_Polycrystal.get_material(1.0,1.0,1.0,init_σ,n0,n1,pixels,grain_size,crystalangle,3,norm_res,norm_mass,alphaN,1.0,mesh,backend)

@testset "JC2D System" begin
    system = TDGL_Polycrystal.jc2d_system(n1,n2,mesh,params,material,backend,B_app)

    @test isa(system.dc,RectPrimalForm2Data2)

    @test system.m.periodic[1] == true
    @test system.m.periodic[2] == yperiodic
    if !yperiodic
        @test sum(system.c.data) == 0
    else
        @test sum(system.c.data) > 0
    end
end
const system = TDGL_Polycrystal.jc2d_system(n1,n2,mesh,params,material,backend,B_app)

@testset "Solver" begin
    solver = TDGL_Polycrystal.ImplicitLondonMultigridSolver(system, state, tol, levelcount)

    @test length(solver.state_u) == levelcount
    @test typeof(solver.rhs[1]) == typeof(state)
end
solver = TDGL_Polycrystal.ImplicitLondonMultigridSolver(system, state, tol, levelcount)

@testset "Jc2D Finder" begin
    finder = TDGL_Polycrystal.JC2DFinder(solver,Ecrit,init_hold_time,wait_time,0.0,Jramp,B_app)

    @test isa(finder,JC2DFinder)
    @test finder.mode == TDGL_Polycrystal.JC2DInitHold()
    @test typeof(finder.esum) == typeof(Ecrit)
    @test finder.solver == solver
    @test sum(finder.δda_rhs.data) == 0
end
finder = TDGL_Polycrystal.JC2DFinder(solver,Ecrit,init_hold_time,wait_time,0.0,Jramp,B_app)

@testset "Boundary Conditions" begin
    f2 = TDGL_Polycrystal.JC2DFinder(solver,Ecrit,init_hold_time,wait_time,0.0,Jramp,B_app)
    
    TDGL_Polycrystal.jc2d_bcs!(finder,system)
    TDGL_Polycrystal.jc2d_bcs!(f2,system)
    TDGL_Polycrystal.jc2d_bcs!(f2,system)

    @test f2.δda_rhs.data == finder.δda_rhs.data

    if !yperiodic
        @test sum(abs.(finder.δda_rhs.data)) > 0
        @test sum(finder.δda_rhs.data) == 0
    else 
        @test sum(abs.(finder.δda_rhs.data)) == 0
    end
end

@testset "Simulation Setup" begin
    FindType = JC2DFinder
    f_jc,metadata,start_α,start_m⁻¹,start_σ = simulation_setup(pixels,N,num_crystal,grain_thick,tstep,GL,init_σ,
    norm_res,norm_mass,Ecrit,Jramp,wait_time,init_hold_time,xmin,ymin,yperiodic,alphaN,betaN,FindType,levelcount,tol,bknd,B_app)

    @test isa(f_jc,FindType)
    @test metadata["Crystal angle"] == atan(1/N)
    @test metadata["Multigrid steps"] == levelcount

    @testset "Step Jc2D" begin
        TDGL_Polycrystal.step!(f_jc)

        @test f_jc.mode == TDGL_Polycrystal.JC2DInitHold()
        @test f_jc.curholdsteps == 1
    end
end

@testset "LinX" begin
    FindType = BVarLinXFinder
    f_jc,metadata,start_α,start_m⁻¹,start_σ = simulation_setup(pixels,N,num_crystal,grain_thick,tstep,GL,init_σ,
    norm_res,norm_mass,Ecrit,Jramp,wait_time,init_hold_time,xmin,ymin,yperiodic,alphaN,betaN,FindType,levelcount,tol,bknd,B_app,Bmax,B_step,num_samples)

    @test isa(f_jc,FindType)
    @test size(f_jc.allJ) == (1+round(abs((B_app-Bmax)/B_step)),num_samples+1)
    @test f_jc.B_index == 1
    @test f_jc.sample_index == 2

    @testset "Step LinX" begin
        TDGL_Polycrystal.step!(f_jc)
        @test f_jc.mode == TDGL_Polycrystal.JC2DInitHold()

        f_jc.mode = TDGL_Polycrystal.BVarLinX()
        f_jc.curholdsteps = 5000
        f_jc.Ehist = range(0.00001,0.0001,f_jc.curholdsteps)
        f_jc.num_holds = f_jc.max_holds
        f_jc.sample_index = f_jc.num_samples + 1
        TDGL_Polycrystal.step!(f_jc)
        @test f_jc.B_field > Bstart
        @test f_jc.mode == TDGL_Polycrystal.BVarLinX()   
    end  
end

@testset "Periodic Boundary Conditions" begin
    backnd = CPU()
    msh = TDGL_Polycrystal.get_mesh(1,xmin,ymin,true)

    c = RectPrimalForm1Data(msh, KernelAbstractions.zeros(backnd, Float64, elemcount(msh,Val(1))))
    dc = RectPrimalForm2Data2(msh, KernelAbstractions.zeros(backnd, Float64, elemcount(msh,Val(2))))

    TDGL_Polycrystal.set_fluxons!(c,dc,0.5,msh,backnd)
    @test sum(asarray(dc, 1, 2)[:,end]) > 0
    @test asarray(c, 2)[1,end] == 0
    @test asarray(c, 2)[end,end] > 0
end


@testset "Simulations" begin 
    @testset "Run Jc2D Finder" begin
        JC2D_finder = TDGL_Polycrystal.JC2DFinder(solver,Ecrit,init_hold_time,wait_time,0.0,Jramp,B_app)
        data, time = TDGL_Polycrystal.find_jc(JC2D_finder,Verbose)

        #@test time < 3
        @test data[2][end] > Ecrit
        @test data[1][end] < 10*Jramp 
        @test data[4][end] == string(TDGL_Polycrystal.JC2DJHold())
    end


    @testset "LinX Reached Equilibrium" begin
        f_linx = BVarLinXFinder(solver,Ecrit,init_hold_time,wait_time,0.0,Jramp,Bstart,Bmax,B_step,num_samples)
        f_linx.mode = TDGL_Polycrystal.BVarLinX()

        @testset "Step J" begin
            #step up J when sample_index = 2
            f_linx.Ehist = range(0.00009,0.0001,1000)
            f_linx.B_index = 1
            TDGL_Polycrystal.reached_equilibrium!(f_linx)
            @test f_linx.sample_index == 3
            @test isapprox(f_linx.allJ[f_linx.B_index,1],Bstart,atol=0.01)
            @test f_linx.B_index == 1
        end

        @testset "Linear Extrapolate" begin
            #linear_extrapolate when sample_index = 3
            f_linx.allJ[1,2] = 0.0001
            f_linx.allE[1,2] = 0.00005
            f_linx.j = 0.0002
            f_linx.Ehist = range(0.000099,0.0001,150)
            TDGL_Polycrystal.reached_equilibrium!(f_linx)
            @test f_linx.j < 0.0002
            @test f_linx.sample_index == 4
        end

        @testset "Out of Bounds" begin
            #nothing happens because an E field val is negative
            f_linx.Ehist = range(-0.00001,-0.0001,1000)
            f_linx.allJ[1,3] = 0.15
            f_linx.allE[1,3] = 0.05
            f_linx.j = 0.2
            f_linx.sample_index = 4
            TDGL_Polycrystal.reached_equilibrium!(f_linx)
            @test f_linx.j == 0.2
            @test f_linx.sample_index == 4
        end
        
        @testset "Next B field" begin
            #move to next B field val when allE/allJ row is full
            f_linx.Ehist = range(0.00001,0.0001,150)
            f_linx.num_holds = 10
            f_linx.sample_index = num_samples+1
            TDGL_Polycrystal.reached_equilibrium!(f_linx)
            @test f_linx.sample_index == 2
            @test f_linx.B_field == f_linx.Brelstep
            @test f_linx.B_index == 2
            @test f_linx.allE[f_linx.B_index,1] == Bstart + f_linx.Brelstep
            @test length(f_linx.Ehist) == 0
            @test f_linx.allJ[f_linx.B_index-1,1] > 0
        end

        @testset "Use Previous B Field" begin
            #extrapolate 2nd j value using gradient from prev B field
            f_linx.Ehist = range(0.00009,0.0001,1000)
            TDGL_Polycrystal.reached_equilibrium!(f_linx)
            @test f_linx.j > 0
        end

        @testset "Simulation Complete" begin
            #finish sim when the full allE array is finished
            f_linx.Ehist = range(0.00009,0.0001,1000)
            f_linx.sample_index = num_samples+1
            f_linx.B_index = f_linx.num_Bs
            TDGL_Polycrystal.reached_equilibrium!(f_linx)
            @test f_linx.mode == TDGL_Polycrystal.JC2DDone()
        end
    end
end


