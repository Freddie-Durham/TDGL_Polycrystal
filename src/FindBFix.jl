
"Run simulation at a fixed J (≈Jc) and B field to study fluctuations in E field"
mutable struct Bfixed{R,VR,VC} <: Finder
    solver::ImplicitLondonMultigridSolver{2,R,RectPrimalForm0Data{2,R,VR},RectPrimalForm1Data{2,R,VR},RectPrimalForm0Data{2,Complex{R},VC},RectPrimalForm1Data{2,Complex{R},VC}}
    mode::LinXMode
    ecrit::R
    shortholdtime::Int64
    longholdtime::Int64
    timesteps::Int64
    curholdsteps::Int64
    max_steps::Int64
    j::R
    jrelstep::R
    E_field::R
    B_field::R
    δda_rhs::RectPrimalForm1Data{2,R,VR}
end

function Bfixed(solver, ecrit::R, shortholdtime, longholdtime, jinit::R, jrelstep::R, startB::R, max_steps) where {R}
    mode = JC2DInitHold()
    timesteps = 0
    curholdsteps = 0
    e_field = zero(R)
    δda_rhs = MulTDGL.similar(MulTDGL.state(solver).a)
    data(δda_rhs) .= zero(eltype(data(δda_rhs)))

    Bfixed(solver,
        mode,
        ecrit,
        shortholdtime,
        longholdtime,
        timesteps,
        curholdsteps,
        max_steps,
        jinit,
        abs(jrelstep),
        e_field,
        startB,
        δda_rhs,
        )
end

function step!(finder::Bfixed)
    sys = system(finder)
    parameters = sys.p

    #update boundary conditions
    jc2d_bcs!(finder,sys)

    #call london multigrid
    step_data = MulTDGL.step!(finder.solver, finder.δda_rhs, (finder.j,0.0)) 
    
    finder.E_field = step_data.e[1]

    finder.timesteps += 1
    finder.curholdsteps += 1

    if finder.mode == JC2DInitHold()
        if finder.curholdsteps >= finder.shortholdtime / parameters.k
            finder.curholdsteps = 0
            finder.mode = JC2DJHold()
        end
    elseif finder.mode == JC2DJHold()
        if finder.E_field < finder.ecrit
            finder.j += finder.jrelstep
            finder.curholdsteps = 0
        elseif finder.curholdsteps >= finder.longholdtime / parameters.k
            finder.mode = BVarLinX()
        end
    elseif finder.mode == BVarLinX() 
        if finder.curholdsteps > finder.max_steps
            finder.mode = JC2DDone()
        end
    end
end

function find_jc(f_jc::Bfixed,verbose::Bool=true)
    current = Vector{Float64}([])
    b_field = Vector{Float64}([])
    e_field = Vector{Float64}([])

    starttime = time()
    while f_jc.mode != JC2DDone()
        push!(b_field,f_jc.B_field)
        push!(current,f_jc.j)
        push!(e_field,f_jc.E_field)
        if verbose
            @time step!(f_jc)
            println("Time taken = $(time()-starttime)")
            println("Current = $(f_jc.j)")
            println("Electric Field = $(f_jc.E_field)")
            println("Magnetic Field = $(f_jc.B_field)")
            println("Current hold: $(f_jc.curholdsteps)")
        else
            step!(f_jc)
        end
    end
    timetaken = time()-starttime

    return [current,b_field,e_field], timetaken
end