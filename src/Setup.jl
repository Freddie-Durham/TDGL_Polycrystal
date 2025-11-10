"applies boundary conditions on every timestep, modifying finder. depends on whether boundaries are periodic in both directions or only one"
function jc_bcs!(finder::Finder,sys::System)
    if sys.m.periodic[1] && !sys.m.periodic[2] && length(sys.m.periodic) == 2
        # for non-periodic systems, a change in the applied field
        # between the upper and lower boundaries results in a net
        # current
        db = finder.j * sys.m.extent[2] * sys.m.h[2] / sys.p.κ^2
        asarray(finder.δda_rhs, 1)[:, 1] .= (finder.B_field - db/2) * 2*sys.m.h[1] / sys.m.h[2]
        asarray(finder.δda_rhs, 1)[:, end] .= -(finder.B_field + db/2) * 2*sys.m.h[1] / sys.m.h[2]

    elseif all(sys.m.periodic)
        # for periodic systems, there can't be a net change in the
        # field. instead, we can apply an external current in the
        # opposite direction so the desired current must flow in the
        # superconductor to compensate
        asarray(finder.δda_rhs, 1) .= -finder.j * sys.m.h[1] / sys.p.κ^2 

    else
        error("Boundary conditions not supported for x non-periodic or not fully periodic 3D systems")
    end
end

"Sets gauge function to fix discontinuity in vector potential due to applied magnetic field in full periodic BCs using fluxon number"
function set_fluxons!(c,dc,B,mesh,backend)
    #flux = fluxoncount * 2π
    #B field = flux / area
    #assume B field always in z direction for now

    flux = B * mesh.extent[1] * mesh.h[1] * mesh.extent[2] * mesh.h[2]

    fluxoncount = round(flux/2π)
    #println("Fluxoncount = $(Int(fluxoncount))")
    n = elemextent(mesh, (2,), 1) #number of y pointing edges in x direction

    if length(mesh.periodic) > 2
        #apply to y pointing edges at top y boundary for all x, z (fixed in z direction, periodic in x direction)
        asarray(c, 2)[:, end, :] .= adapt(backend, [fluxoncount * (i - 1) / n  * 2π 
        for i in 1:n, _ in 1:elemextent(mesh, (2,), 3)])
        asarray(dc, 1, 2)[:, end, :] .= fluxoncount / n * 2π

    else #for 2D periodic systems
        asarray(c, 2)[:, end] .= adapt(backend, [fluxoncount * (i - 1) / n  * 2π for i in 1:n])
        asarray(dc, 1, 2)[:, end] .= fluxoncount / n * 2π

    end
end

"Creates system with phase difference and parallel transport. Different for doubly periodic BCs"
function jc_system(n1, n2, mesh::RectMesh{N,R}, params, material, backend, B_field) where {N,R}
    # phase difference due to change of local section
    c = RectPrimalForm1Data(mesh, KernelAbstractions.zeros(backend, R, n1))
   
    if length(mesh.periodic) > 2 # 3D system
        dc = RectPrimalForm2Data3(mesh, KernelAbstractions.zeros(backend, R, n2))
    else # 2D system
        dc = RectPrimalForm2Data2(mesh, KernelAbstractions.zeros(backend, R, n2))
    end

    if all(mesh.periodic)
        # the wavefunction is a section of a fibre bundle and the
        # vector potential is essentially the connection. by putting a
        # twist in the bundle we can have any whole number of flux
        # quanta even in a periodic domain (a flat torus)
        set_fluxons!(c,dc,B_field,mesh,backend)
    end

    # parallel transport (updated automatically)
    u = RectPrimalForm1Data(mesh, KernelAbstractions.ones(backend, Complex{R}, n1))

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

"Returns a state for a system"
function get_state(n0,n1,m::RectMesh{N,R},backend,seed) where {N,R}
    # initial state
    ψdata = Adapt.adapt(backend, exp.(Complex{R}(2π *im) .* rand(Xoshiro(seed), R, n0)))
    ψ = RectPrimalForm0Data(m, ψdata)        # wavefunction with random phase (same every time function is called)
    a = RectPrimalForm1Data(m, KernelAbstractions.zeros(backend, R, n1))
    φ = RectPrimalForm0Data(m, KernelAbstractions.zeros(backend, R, n0))

    return MulTDGL.State(ψ, a, φ)
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


"Returns a rectilinear mesh for a 2D system, accounting for requirements of multigrid"
function get_mesh(pixels,xdim,ydim,yperiodic)
    # simulation mesh
    extent = (xdim,ydim).*pixels
    periodic = (true, yperiodic)
    h = (1.0/pixels, 1.0/pixels) # grid step
    return MulTDGL.RectMesh(extent, periodic, h)
end

"Returns a rectilinear mesh for a 3D system, accounting for requirements of multigrid"
function get_3D_mesh(type, pixels, xdim, ydim, zdim, yperiodic=true, zperiodic=true)
    # simulation mesh
    extent = (xdim, ydim, zdim) .* pixels
    periodic = (true, yperiodic, zperiodic)
    h = (type(1.0/pixels), type(1.0/pixels), type(1.0/pixels)) # grid step
    return MulTDGL.RectMesh(extent, periodic, h)
end

"Setup material using Crystal2D according to Carty(2008)"
function get_material(init_α,init_β,init_m⁻¹,init_σ,pixels,pattern,norm_resist,norm_mass,alphaN,betaN,m,backend)
    #relative weights of nodes based on whether inside or outside grain boundaries
    weights = apply_pattern(pattern,pixels,elemextent(m, ()))
    #apply pattern to normalised α value
    start_α = map(w->linear_interp(init_α,alphaN,w), weights)

    #apply coating if non-periodic in y direction
    if m.periodic[1] && !m.periodic[2] && typeof(pattern) == TDGL_Polycrystal.TruncOct
        thickness = Int(7.5*pixels)
        start_α[:,1:thickness] .= alphaN
        start_α[:,end-thickness:end] .= alphaN
    end
    α = RectPrimalForm0Data(m, adapt(backend, reshape(start_α, :)))

    start_β = map(w->linear_interp(init_β,betaN,w), weights)
    β = RectPrimalForm0Data(m, adapt(backend, reshape(start_β, :)))
    
    #need two lots of conductivity for x and y 
    start_σX = map(w->linear_interp(init_σ,1/norm_resist,w),interp_edgesX(weights,m.periodic[1]))
    start_σY = map(w->linear_interp(init_σ,1/norm_resist,w),interp_edgesY(weights,m.periodic[2]))
    σ = RectPrimalForm1Data(m, adapt(backend, vcat(reshape(start_σX, :),reshape(start_σY, :)))) 

    #same goes for reciprocal mass
    start_m⁻¹X = map(w->linear_interp(init_m⁻¹,1/norm_mass,w),interp_edgesX(weights,m.periodic[1]))
    start_m⁻¹Y = map(w->linear_interp(init_m⁻¹,1/norm_mass,w),interp_edgesY(weights,m.periodic[2]))
    m⁻¹ = RectPrimalForm1Data(m, adapt(backend, vcat(reshape(start_m⁻¹X, :),reshape(start_m⁻¹Y, :)))) 

    return MulTDGL.Material(α, β, m⁻¹, σ),start_α,start_β,start_m⁻¹X,start_σY
end

"Setup material using Voronoi according to Blair(2022)"
function get_3D_material(init_α,init_β,init_m⁻¹,init_σ,pixels,pattern,norm_resist,norm_mass,alphaN,betaN,m::RectMesh{N,R},backend) where {N,R}
    #relative weights of nodes based on whether inside or outside grain boundaries
    weights = apply_3D_pattern(pattern,pixels,elemextent(m, ()))

    #apply pattern to normalised α value
    start_α = map(w->R(linear_interp(init_α,alphaN,w)), weights)
    α = RectPrimalForm0Data(m, adapt(backend, reshape(start_α, :)))

    start_β = map(w->R(linear_interp(init_β,betaN,w)), weights)
    β = RectPrimalForm0Data(m, adapt(backend, reshape(start_β, :)))

    weights_x,weights_y,weights_z = interpolate_edges(weights,m.periodic)

    σ = RectPrimalForm1Data(m, adapt(backend, vcat(
        reshape(map(w->R(linear_interp(init_σ,1/norm_resist,w)),weights_x),:),
        reshape(map(w->R(linear_interp(init_σ,1/norm_resist,w)),weights_y),:),
        reshape(map(w->R(linear_interp(init_σ,1/norm_resist,w)),weights_z),:)
    )))

    m⁻¹ = RectPrimalForm1Data(m, adapt(backend, vcat(
        reshape(map(w->R(linear_interp(init_m⁻¹,1/norm_mass,w)),weights_x),:),
        reshape(map(w->R(linear_interp(init_m⁻¹,1/norm_mass,w)),weights_y),:),
        reshape(map(w->R(linear_interp(init_m⁻¹,1/norm_mass,w)),weights_z),:)
    )))

    return MulTDGL.Material(α, β, m⁻¹, σ), weights
end

"Set up parameters of simulation using CUDA"
function simulation_setup(vortex_radius,pattern,tstep,GL,init_σ,norm_resist,norm_mass,Ecrit,Jramp,wait_time,init_hold_time,xdim,ydim,yperiodic,alphaN,betaN,init_alpha,init_beta,finder,levelcount,tol,bknd,rng_seed,J_init,B_init,args...)
    backend = get_backend(bknd)

    #get parameters + material, setup the SC system. Create initial state of system and initialise solver. Return solver along with BCs
    params = get_params(tstep,GL)
    mesh = get_mesh(vortex_radius,xdim,ydim,yperiodic)

    # number of elements of given degree
    n0 = elemcount(mesh, Val(0)) # number of vertices
    n1 = elemcount(mesh, Val(1)) # number of edges
    n2 = elemcount(mesh, Val(2)) # number of faces

    init_m⁻¹ = 1.0
    material,start_α,start_β,start_m⁻¹X,start_σY = get_material(
        init_alpha,init_beta,init_m⁻¹,init_σ,vortex_radius,pattern,norm_resist,norm_mass,alphaN,betaN,mesh,backend)
    
    system = jc_system(n1,n2,mesh,params,material,backend,B_init)
    state = get_state(n0,n1,mesh,backend,rng_seed)

    # create simulation
    s = MulTDGL.ImplicitLondonMultigridSolver(system, state, tol, levelcount)

    f_jc = finder(s,Ecrit,init_hold_time,wait_time,J_init,Jramp,B_init,args...)

    metadata = Dict("Timestep" => tstep, "κ" => GL, "xMesh" => mesh.extent[1],
    "yMesh" => mesh.extent[2], "Vortex radius (pixels)" => vortex_radius, "xPeriodic" => string(mesh.periodic[1]),
    "yPeriodic" => string(mesh.periodic[2]), "Tolerance" => tol, "Multigrid steps" => levelcount,
    "α_s" => init_alpha, "β_s" => init_beta, "Effective electron mass" => norm_mass, "Normal resistivity" => norm_resist,
    "α_n" => alphaN, "β_n" => betaN,"Initial hold time" => init_hold_time, "Wait to stabilise" => wait_time,
    "J ramp" => Jramp, "E criterion" => Ecrit,
    "Date" => string(now()), "Backend" => bknd, "Random seed" => rng_seed,"Finder" => string(finder),"MulTDGL Version no." => string(pkgversion(MulTDGL)),
    "TDGL_Polycrystal Version no." => Version)
    
    append_metadata!(metadata,pattern)

    return f_jc,metadata,start_α,start_β,start_m⁻¹X,start_σY
end

"Set up parameters of simulation using CUDA"
function simulation_setup_3D(vortex_radius,pattern,tstep,GL,init_σ,norm_resist,norm_mass,Ecrit,Jramp,wait_time,init_hold_time,xdim,ydim,zdim,yperiodic,zperiodic,alphaN,betaN,init_alpha,init_beta,finder,levelcount,tol,bknd,rng_seed,J_init,B_init,args...)
    type = Float32
    backend = get_backend(bknd)
    params = get_params(type(tstep),type(GL))
    mesh = get_3D_mesh(type,vortex_radius,xdim,ydim,zdim,yperiodic,zperiodic)

    # number of elements of given degree
    n0 = elemcount(mesh, Val(0)) # number of vertices
    n1 = elemcount(mesh, Val(1)) # number of edges
    n2 = elemcount(mesh, Val(2)) # number of faces

    init_m⁻¹ = 1.0
    material,weights = get_3D_material(
    init_alpha,init_beta,init_m⁻¹,init_σ,vortex_radius,pattern,
    norm_resist,norm_mass,alphaN,betaN,mesh,backend)

    system = jc_system(n1,n2,mesh,params,material,backend,B_init)
    state = get_state(n0,n1,mesh,backend,rng_seed)

    s = MulTDGL.ImplicitLondonMultigridSolver(system, state, type(tol), levelcount)
    f_jc = finder(s,Ecrit,init_hold_time,wait_time,J_init,Jramp,B_init,args...)

    metadata = Dict("Timestep" => tstep, "κ" => GL, "xMesh" => mesh.extent[1],
    "yMesh" => mesh.extent[2], "zMesh" => mesh.extent[3], "Vortex radius (pixels)" => vortex_radius, 
    "Periodic" => string(mesh.periodic), "Tolerance" => tol, "Multigrid steps" => levelcount,
    "α_s" => init_alpha, "β_s" => init_beta, "α_n" => alphaN, "β_n" => betaN,
    "GB Effective Mass" => norm_mass, "GB Resistivity" => norm_resist, "SC Effective Mass" => 1.0/init_m⁻¹, "SC Resistivity" => 1.0/init_σ,
    "Initial hold time" => init_hold_time, "Wait to stabilise" => wait_time,
    "J ramp" => Jramp, "E criterion" => Ecrit,"Date" => string(now()),
    "Backend" => bknd, "Random seed" => rng_seed,"Finder" => string(finder),
    "MulTDGL Version no." => string(pkgversion(MulTDGL)), "TDGL_Polycrystal Version no." => Version)
    
    append_metadata!(metadata,pattern)

    return f_jc,metadata,weights
end

"set the B field of the system in the solver and reset the solver state"
function new_solver(solver::ImplicitLondonMultigridSolver{N,R,VR,VC},B_field,tol,levelcount,backend,rng) where {N,R,VR,VC}
    old_sys = MulTDGL.system(solver)
    mesh = old_sys.m
    n0 = elemcount(mesh, Val(0)) # number of vertices
    n1 = elemcount(mesh, Val(1)) # number of edges
    n2 = elemcount(mesh, Val(2)) # number of faces
    state = get_state(n0,n1,mesh,backend,rng)
    system = jc_system(n1,n2,mesh,old_sys.p,old_sys.mat,backend,B_field)
    return MulTDGL.ImplicitLondonMultigridSolver(system, state, R(tol), levelcount)
end

"Create a new finder for a new B field and applied current density"
function new_finder(F::Finder,findtype,Ecrit,wait_time,init_hold_time,Jramp,J_init,B_init,tol,levelcount,backend,rng,args...)
    s = new_solver(F.solver,B_init,tol,levelcount,get_backend(backend),rng)
    return findtype(s,Ecrit,init_hold_time,wait_time,J_init,Jramp,B_init,args...)
end

