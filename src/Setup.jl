"applies boundary conditions on every timestep, modifying finder. depends on whether boundaries are periodic in both directions or only one"
function jc2d_bcs!(finder::Finder,sys::System)
    if sys.m.periodic[1] && !sys.m.periodic[2]
        # for non-periodic systems, a change in the applied field
        # between the upper and lower boundaries results in a net
        # current
        db = finder.j * sys.m.extent[2] * sys.m.h[2] / sys.p.κ^2
        asarray(finder.δda_rhs, 1)[:, 1] .= (finder.B_field - db/2) * 2*sys.m.h[1] / sys.m.h[2]
        asarray(finder.δda_rhs, 1)[:, end] .= -(finder.B_field + db/2) * 2*sys.m.h[1] / sys.m.h[2]

    elseif sys.m.periodic[1] && sys.m.periodic[2]
        # for periodic systems, there can't be a net change in the
        # field. instead, we can apply an external current in the
        # opposite direction so the desired current must flow in the
        # superconductor to compensate
        asarray(finder.δda_rhs, 1) .= -finder.j * sys.m.h[1] / sys.p.κ^2 
    end
end

"Sets gauge function to fix discontinuity in vector potential due to applied magnetic field in full periodic BCs using fluxon number"
function set_fluxons!(c,dc,B,mesh,backend)
    #flux = fluxoncount * 2π
    #B field = flux / area
    flux = B * mesh.extent[1] * mesh.h[1] * mesh.extent[2] * mesh.h[2]

    fluxoncount = round(flux/2π)
    println("Fluxoncount = $(Int(fluxoncount))")

    n = elemextent(mesh, (2,), 1)
    asarray(c, 2)[:, end] .=
        adapt(backend, [fluxoncount * (i - 1) / n  * 2π for i in 1:n])

    asarray(dc, 1, 2)[:, end] .= fluxoncount / n * 2π
end

"Creates system with phase difference and parallel transport. Different for doubly periodic BCs"
function jc2d_system(n1, n2, mesh, params, material, backend, B_field)
    # phase difference due to change of local section
    c = RectPrimalForm1Data(mesh, KernelAbstractions.zeros(backend, Float64, n1))
    dc = RectPrimalForm2Data2(mesh, KernelAbstractions.zeros(backend, Float64, n2))

    if mesh.periodic[1] && mesh.periodic[2]
        #println("Doubly periodic BCs")
        # the wavefunction is a section of a fibre bundle and the
        # vector potential is essentially the connection. by putting a
        # twist in the bundle we can have any whole number of flux
        # quanta even in a periodic domain (a flat torus)
        set_fluxons!(c,dc,B_field,mesh,backend)
    end

    # parallel transport (updated automatically)
    u = RectPrimalForm1Data(mesh, KernelAbstractions.ones(backend, ComplexF64, n1))

    # create system
    return System(mesh, params, material, c, dc, u)
end


"returns parameters timestep (<1) and kappa (>1 for technological SCs)"
function get_params(tstep,GL)
    # parameters
    k = tstep # timestep                 
    κ = GL    # Ginzburg–Landau parameter
   return MulTDGL.Parameters(k, κ)
end


"Returns a rectilinear mesh for a 2D system, accounting for requirements of multigrid"
function get_mesh(pixels,xdim,ydim,yperiodic)
    # simulation mesh
    extent = (xdim,ydim).*pixels
    periodic = (true, yperiodic)
    h = (1.0/pixels, 1.0/pixels) # grid step
    return MulTDGL.RectMesh(extent, periodic, h)
end

"Returns a state for a 2D system"
function get_state(n0,n1,m,backend,seed)
    # initial state
    ψdata = Adapt.adapt(backend, exp.(2π * im .* rand(Xoshiro(seed), Float64, n0)))
    ψ = RectPrimalForm0Data(m, ψdata)        # wavefunction with random phase (same every time function is called)
    a = RectPrimalForm1Data(m, KernelAbstractions.zeros(backend, Float64, n1))
    φ = RectPrimalForm0Data(m, KernelAbstractions.zeros(backend, Float64, n0))

    return MulTDGL.State(ψ, a, φ)
end

"Setup material using Crystal2D according to Carty(2008) - 2D slice of truncated octahedra"
function get_material(init_α,init_β,init_m⁻¹,init_σ,pixels,factor,grain_size,crystalangle,grain_thick,norm_resist,norm_mass,alphaN,betaN,m,backend)
    #apply pattern to normalised α value
    start_α = apply_pattern(octagon!,init_α,alphaN,pixels,factor,grain_size,grain_thick,crystalangle,m.extent[1],m.extent[2],elemextent(m, ()))
    #apply coating if non-periodic in y direction
    if m.periodic[1] && !m.periodic[2]
        thickness = Int(7.5*pixels)
        start_α[:,1:thickness] .= alphaN
        start_α[:,end-thickness:end] .= alphaN
    end
    α = RectPrimalForm0Data(m, adapt(backend, reshape(start_α, :)))
    
    #need two lots of conductivity for x and y 
    start_σ = apply_pattern(octagon!,init_σ,1/norm_resist,pixels,factor,grain_size,grain_thick,crystalangle,m.extent[1],m.extent[2],elemextent(m, (1,)))
    start_2σ = apply_pattern(octagon!,init_σ,1/norm_resist,pixels,factor,grain_size,grain_thick,crystalangle,m.extent[1],m.extent[2],elemextent(m, (2,)))
    σ = RectPrimalForm1Data(m, adapt(backend, vcat(reshape(start_σ, :),reshape(start_2σ, :)))) 

    #same goes for reciprocal mass
    start_m⁻¹ = apply_pattern(octagon!,init_m⁻¹,1/norm_mass,pixels,factor,grain_size,grain_thick,crystalangle,m.extent[1],m.extent[2],elemextent(m, (1,)))
    start_2m⁻¹ = apply_pattern(octagon!,init_m⁻¹,1/norm_mass,pixels,factor,grain_size,grain_thick,crystalangle,m.extent[1],m.extent[2],elemextent(m, (2,)))
    m⁻¹ = RectPrimalForm1Data(m, adapt(backend, vcat(reshape(start_m⁻¹, :),reshape(start_2m⁻¹, :)))) 

    #can also change β but often not necessary
    start_β = apply_pattern(octagon!,init_β,betaN,pixels,factor,grain_size,grain_thick,crystalangle,m.extent[1],m.extent[2],elemextent(m, ()))
    β = RectPrimalForm0Data(m, adapt(backend, reshape(start_β, :)))
    
    return MulTDGL.Material(α, β, m⁻¹, σ),start_α,start_β,start_m⁻¹,start_σ
end

"Set up single Josephson junction material pattern"
function set_material(init_α,init_β,init_m⁻¹,init_σ,pixels,factor,grain_size,crystalangle,junc_thick,norm_resist,norm_mass,alphaN,betaN,m,backend)
    #apply pattern to normalised α value
    start_α = apply_pattern(init_α,alphaN,pixels,factor,junc_thick,m.extent[1],m.extent[2],elemextent(m, ()))
    α = RectPrimalForm0Data(m, adapt(backend, reshape(start_α, :)))
    
    #need two lots of conductivity for x and y 
    start_σ = apply_pattern(init_σ,1/norm_resist,pixels,factor,junc_thick,m.extent[1],m.extent[2],elemextent(m, (1,)))
    start_2σ = apply_pattern(init_σ,1/norm_resist,pixels,factor,junc_thick,m.extent[1],m.extent[2],elemextent(m, (2,)))
    σ = RectPrimalForm1Data(m, adapt(backend, vcat(reshape(start_σ, :),reshape(start_2σ, :)))) 

    #same goes for reciprocal mass
    start_m⁻¹ = apply_pattern(init_m⁻¹,1/norm_mass,pixels,factor,junc_thick,m.extent[1],m.extent[2],elemextent(m, (1,)))
    start_2m⁻¹ = apply_pattern(init_m⁻¹,1/norm_mass,pixels,factor,junc_thick,m.extent[1],m.extent[2],elemextent(m, (2,)))
    m⁻¹ = RectPrimalForm1Data(m, adapt(backend, vcat(reshape(start_m⁻¹, :),reshape(start_2m⁻¹, :)))) 

    #can also change β but often not necessary
    start_β = apply_pattern(init_β,betaN,pixels,factor,junc_thick,m.extent[1],m.extent[2],elemextent(m, ()))
    β = RectPrimalForm0Data(m, adapt(backend, reshape(start_β, :)))
    
    return MulTDGL.Material(α, β, m⁻¹, σ),start_α,start_β,start_m⁻¹,start_σ
end

"Convert string to type of backend"
function get_backend(bknd)
    if uppercase(bknd) == "CUDA"
        return CUDABackend() # for NVIDIA GPUs
    elseif uppercase(bknd) == "AMDGPU"
        return AMDGPUBackend() # for AMD GPUs
    elseif uppercase(bknd) == "CPU"
        return CPU() # for CPU
    else
        error("Backend:$bknd not found")
    end
end

"Set up parameters of simulation using CUDA"
function simulation_setup(vortex_radius,factor,N,rep_grain,grain_thick,tstep,GL,init_σ,norm_resist,norm_mass,Ecrit,Jramp,wait_time,init_hold_time,xdim,ydim,yperiodic,alphaN,betaN,init_alpha,init_beta,finder,levelcount,tol,bknd,rng_seed,B_init,args...)
    backend = get_backend(bknd)

    #setup grain size and orientation to ensure periodicity
    grainangle,grain_diameter = periodic_crystal(N,fld(xdim,rep_grain))

    #get parameters + material, setup the SC system. Create initial state of system and initialise solver. Return solver along with BCs
    params = get_params(tstep,GL)
    mesh = get_mesh(vortex_radius,xdim,ydim,yperiodic)

    # number of elements of given degree
    n0 = elemcount(mesh, Val(0)) # number of vertices
    n1 = elemcount(mesh, Val(1)) # number of edges
    n2 = elemcount(mesh, Val(2)) # number of faces

    init_m⁻¹ = 1.0
    material,start_α,start_β,start_m⁻¹,start_σ = get_material(init_alpha,init_beta,init_m⁻¹,init_σ,vortex_radius,factor,grain_diameter,grainangle,grain_thick,norm_resist,norm_mass,alphaN,betaN,mesh,backend)
    #material,start_α,start_β,start_m⁻¹,start_σ = set_material(init_alpha,init_beta,init_m⁻¹,init_σ,vortex_radius,factor,grain_diameter,grainangle,grain_thick,norm_resist,norm_mass,alphaN,betaN,mesh,backend)

    system = jc2d_system(n1,n2,mesh,params,material,backend,B_init)
    state = get_state(n0,n1,mesh,backend,rng_seed)

    # create simulation
    s = MulTDGL.ImplicitLondonMultigridSolver(system, state, tol, levelcount)

    f_jc = finder(s,Ecrit,init_hold_time,wait_time,0.0,Jramp,B_init,args...)

    metadata = Dict("Timestep" => tstep, "κ" => GL, "xMesh" => mesh.extent[1],
    "yMesh" => mesh.extent[2], "Vortex radius (pixels)" => vortex_radius, "xPeriodic" => string(mesh.periodic[1]),
    "yPeriodic" => string(mesh.periodic[2]), "Tolerance" => tol, "Multigrid steps" => levelcount,
    "α_s" => init_alpha, "β_s" => init_beta, "Effective electron mass" => norm_mass, "Normal resistivity" => norm_resist,
    "α_n" => alphaN, "β_n" => betaN,
    "Grain size (ξ)" => grain_diameter,"Multiple of grains" => rep_grain,"Initial hold time" => init_hold_time, "Wait to stabilise" => wait_time,
    "Lattice angle" => grainangle,  "J ramp" => Jramp, "E criterion" => Ecrit, "Grain Boundary Thickness (ξ)" => grain_thick,
    "Date" => string(now()), "Backend" => bknd, "Finder" => string(finder),"MulTDGL Version no." => string(pkgversion(MulTDGL)),
    "TDGL_Polycrystal Version no." => Version)

    return f_jc,metadata,start_α,start_β,start_m⁻¹,start_σ
end

"set the B field of the system in the solver and reset the solver state"
function new_solver(solver,B_field,tol,levelcount,backend,rng)
    old_sys = MulTDGL.system(solver)
    mesh = old_sys.m
    n0 = elemcount(mesh, Val(0)) # number of vertices
    n1 = elemcount(mesh, Val(1)) # number of edges
    n2 = elemcount(mesh, Val(2)) # number of faces
    state = get_state(n0,n1,mesh,backend,rng)
    system = jc2d_system(n1,n2,mesh,old_sys.p,old_sys.mat,backend,B_field)
    return MulTDGL.ImplicitLondonMultigridSolver(system, state, tol, levelcount)
end

"Create a new finder for a new B field and applied current density"
function new_finder(F::Finder,finder,Ecrit,wait_time,init_hold_time,Jramp,B_init,tol,levelcount,backend,rng,args...)
    s = new_solver(F.solver,B_init,tol,levelcount,get_backend(backend),rng)
    return finder(s,Ecrit,init_hold_time,wait_time,0.0,Jramp,B_init,args...)
end

