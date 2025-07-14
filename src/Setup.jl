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

"Sets phase difference to mimic applied magnetic field in full periodic BCs using fluxons"
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
        println("Doubly periodic BCs")
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
    k = tstep # timestep                  0.5
    κ = GL    # Ginzburg–Landau parameter 10.0
   return MulTDGL_FD.Parameters(k, κ)
end


"Returns a rectilinear mesh for a 2D system, accounting for requirements of multigrid"
function get_mesh(pixels,xdim,ydim,yperiodic)
    # simulation mesh
    extent = (xdim,ydim).*pixels
    periodic = (true, yperiodic)
    h = (1.0/pixels, 1.0/pixels) # grid step
    return MulTDGL_FD.RectMesh(extent, periodic, h)
end


"Returns a state for a 2D system"
function get_state(n0,n1,m,backend)
    # initial state
    ψdata = Adapt.adapt(backend, exp.(2π * im .* rand(Float64, n0)))
    ψ = RectPrimalForm0Data(m, ψdata)                          # wavefunction with random phase
    a = RectPrimalForm1Data(m, KernelAbstractions.zeros(backend, Float64, n1))
    φ = RectPrimalForm0Data(m, KernelAbstractions.zeros(backend, Float64, n0))

    return MulTDGL_FD.State(ψ, a, φ)
end


"Setup material using CrystalLattice code according to Carty(2008)"
function get_material(init_α,init_β,init_m⁻¹,init_σ,n0,n1,pixels,grain_size,crystalangle,grain_thick,norm_resist,norm_mass,alphaN,betaN,m,backend)
    start_α = init_α*ones(Float64, elemextent(m, ()))

    if !(m.periodic[1] && m.periodic[2]) #we dont need alpha gradient to reduce edge effects if doubly periodic BCs
        CrystalLattice.non_periodic!(start_α,elemextent(m, ()),1,10*pixels,alphaN,1.0)
    end
    
    CrystalLattice.tesselate!(CrystalLattice.octagon!,start_α,grain_size,grain_thick*pixels,crystalangle,m.extent[1],m.extent[2],alphaN)
    α = RectPrimalForm0Data(m, adapt(backend, reshape(start_α, :)))
    
    #need two lots of resistivity for x and y 
    start_σ = init_σ*ones(Float64, elemextent(m, (1,)))
    start_2σ = init_σ*ones(Float64, elemextent(m, (2,)))
    CrystalLattice.tesselate!(CrystalLattice.octagon!,start_σ,grain_size,grain_thick*pixels,crystalangle,m.extent[1],m.extent[2],1/norm_resist)
    CrystalLattice.tesselate!(CrystalLattice.octagon!,start_2σ,grain_size,grain_thick*pixels,crystalangle,m.extent[1],m.extent[2],1/norm_resist)

    start_m⁻¹ = init_m⁻¹*ones(Float64, elemextent(m, (1,)))
    start_2m⁻¹ = init_m⁻¹*ones(Float64, elemextent(m, (2,)))
    CrystalLattice.tesselate!(CrystalLattice.octagon!,start_m⁻¹,grain_size,grain_thick*pixels,crystalangle,m.extent[1],m.extent[2],1/norm_mass)
    CrystalLattice.tesselate!(CrystalLattice.octagon!,start_2m⁻¹,grain_size,grain_thick*pixels,crystalangle,m.extent[1],m.extent[2],1/norm_mass)

    σ = RectPrimalForm1Data(m, adapt(backend, vcat(reshape(start_σ, :),reshape(start_2σ, :))))       #normal state resistivity
    m⁻¹ = RectPrimalForm1Data(m, adapt(backend, vcat(reshape(start_m⁻¹, :),reshape(start_2m⁻¹, :)))) #reciprocal normal effective mass

    start_β = init_β*ones(Float64, elemextent(m, ()))
    CrystalLattice.tesselate!(CrystalLattice.octagon!,start_β,grain_size,grain_thick*pixels,crystalangle,m.extent[1],m.extent[2],betaN)
    β = RectPrimalForm0Data(m, adapt(backend, reshape(start_β, :)))
    return MulTDGL_FD.Material(α, β, m⁻¹, σ),start_α,start_β,start_m⁻¹,start_σ
end


"Set up parameters of simulation using CUDA"
function simulation_setup(vortex_radius,N,num_crystal,grain_thick,tstep,GL,init_σ,norm_resist,norm_mass,Ecrit,Jramp,wait_time,init_hold_time,xdim,ydim,yperiodic,alphaN,betaN,finder,levelcount,tol,bknd,Version,B_init,args...)
    if bknd == "CUDA"
        backend = CUDABackend() # for NVIDIA GPUs
    elseif bknd == "CPU"
        backend = CPU() # for CPU
    else
        error("Backend:$bknd not found")
    end

    #setup grain size and orientation to ensure periodicity
    crystalangle,crystal_diameter = periodic_crystal(N,xdim)
    crystal_diameter *= 2.0^(-num_crystal)
    grain_size = crystal_diameter*vortex_radius #width in pixels

    #get parameters + material, setup the SC system. Create initial state of system and initialise solver. Return solver along with BCs
    params = get_params(tstep,GL)
    mesh = get_mesh(vortex_radius,xdim,ydim,yperiodic)

    # number of elements of given degree
    n0 = elemcount(mesh, Val(0)) # number of vertices
    n1 = elemcount(mesh, Val(1)) # number of edges
    n2 = elemcount(mesh, Val(2)) # number of faces

    init_α = 1.0
    init_β = 1.0
    init_m⁻¹ = 1.0
    material,start_α,start_β,start_m⁻¹,start_σ = get_material(init_α,init_β,init_m⁻¹,init_σ,n0,n1,vortex_radius,grain_size,crystalangle,grain_thick,norm_resist,norm_mass,alphaN,betaN,mesh,backend)

    system = jc2d_system(n1,n2,mesh,params,material,backend,B_init)
    state = get_state(n0,n1,mesh,backend)

    # create simulation
    s = ImplicitLondonMultigridSolver(system, state, tol, levelcount)

    f_jc = finder(s,Ecrit,init_hold_time,wait_time,0.0,Jramp,B_init,args...)

    metadata = Dict("Timestep" => tstep, "κ" => GL, "xMesh" => mesh.extent[1],
    "yMesh" => mesh.extent[2], "Vortex radius (pixels)" => vortex_radius, "xPeriodic" => string(mesh.periodic[1]),
    "yPeriodic" => string(mesh.periodic[2]), "Tolerance" => tol, "Multigrid steps" => levelcount, 
    "α" => init_α, "β" => init_β, "Effective electron mass" => norm_mass, "Normal resistivity" => norm_resist,
    "Grain size" => grain_size,"Multiple of grains" => 2^num_crystal,"Initial hold time" => init_hold_time, "Wait to stabilise" => wait_time,
    "Crystal angle" => crystalangle,  "J ramp" => Jramp, "E criterion" => Ecrit, "Grain Boundary Thickness" => grain_thick,
    "Date" => string(now()), "Backend" => string(backend), "Finder" => string(finder),"MulTDGL Version no." => string(pkgversion(MulTDGL_FD)),
    "TDGL2D Version no." => string(Version), "CrystalLattice Version no."=> string(pkgversion(CrystalLattice)))

    return f_jc,metadata,start_α,start_β,start_m⁻¹,start_σ
end
