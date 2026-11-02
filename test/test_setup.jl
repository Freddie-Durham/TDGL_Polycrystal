@testset "Parameters & Backend" begin
    params = TDGL_Polycrystal.get_params(tstep, GL)
    @test isa(params, Parameters)
    @test typeof(params.k) == Float64
    @test params.k == tstep

    @test TDGL_Polycrystal.interpret_backend("CPU") isa CPU
    @test TDGL_Polycrystal.interpret_backend("cuda") isa CUDABackend
    @test_throws ErrorException TDGL_Polycrystal.interpret_backend("UNKNOWN_BACKEND")
end

@testset "Mesh & State (2D & 3D)" begin
    mesh2d = TDGL_Polycrystal.get_mesh(pixels, xmin, ymin, false)
    @test isa(mesh2d, RectMesh)
    @test mesh2d.extent == (xmin, ymin) .* pixels
    @test mesh2d.periodic == (true, false)

    n0 = elemcount(mesh2d, Val(0))
    n1 = elemcount(mesh2d, Val(1))
    state2d = TDGL_Polycrystal.get_state(n0, n1, mesh2d, CPU(), 42)
    @test isa(state2d, State)
    @test sum(Adapt.adapt(CPU(), data(state2d.a))) == 0
    @test TDGL_Polycrystal.MulTDGL_backend(state2d.ψ) isa CPU

    mesh3d = TDGL_Polycrystal.get_3D_mesh(Float64, pixels, xmin, ymin, zmin, true, false)
    @test isa(mesh3d, RectMesh)
    @test mesh3d.extent == (xmin, ymin, zmin) .* pixels
    @test mesh3d.periodic == (true, true, false)
end

@testset "Material & System (2D & 3D)" begin
    # 2D JosephsonJunction (non-periodic in y)
    jj_pattern = TDGL_Polycrystal.JosephsonJunction(2.0, 6)
    np_mesh = TDGL_Polycrystal.get_mesh(1, 8, 16, false)
    mat_jj, start_α, start_β, start_m⁻¹X, start_m⁻¹Y = TDGL_Polycrystal.get_material(
        1.0, 1.0, 1.0, init_σ, 1, jj_pattern, norm_res, norm_mass, alphaN, betaN, np_mesh, CPU())
    @test isa(mat_jj, Material)
    @test start_α[1, 1] == alphaN
    @test start_α[1, 8] == 1.0

    # 2D Voronoi (periodic in y) using D500AR100xyID1.tess (grain_size=8.0, pixels=1, dims=(16,16))
    vor2d = TDGL_Polycrystal.Voronoi(8.0, 1.0, 1)
    p_mesh = TDGL_Polycrystal.get_mesh(1, 16, 16, true)
    mat_vor, α_v, β_v, m_x, m_y = TDGL_Polycrystal.get_material(
        1.0, 1.0, 1.0, init_σ, 1, vor2d, norm_res, norm_mass, alphaN, betaN, p_mesh, CPU())
    @test isa(mat_vor, Material)

    params = TDGL_Polycrystal.get_params(tstep, GL)
    sys_np = TDGL_Polycrystal.jc_system(elemcount(np_mesh, Val(1)), elemcount(np_mesh, Val(2)), np_mesh, params, mat_jj, CPU(), B_app)
    @test sum(Adapt.adapt(CPU(), data(sys_np.c))) == 0

    sys_p = TDGL_Polycrystal.jc_system(elemcount(p_mesh, Val(1)), elemcount(p_mesh, Val(2)), p_mesh, params, mat_vor, CPU(), B_app)
    @test sum(Adapt.adapt(CPU(), data(sys_p.c))) > 0
end

@testset "Simulation Setup & Boundary Conditions (2D & 3D)" begin
    vor2d = TDGL_Polycrystal.Voronoi(8.0, 1.0, 1)
    # 2D Periodic setup
    f_2d_p, meta_2d_p, s_α, s_β, s_mx, s_my = simulation_setup(
        1, vor2d, tstep, GL, init_σ, norm_res, norm_mass, Ecrit, Jramp,
        wait_time, init_hold_time, 16, 16, true, alphaN, betaN, 1.0, 1.0,
        JcFinder, 2, tol, "CPU", 1, 1e-4, B_app, "EXP_INCREASE")
    @test f_2d_p isa JcFinder
    @test meta_2d_p["Pattern"] == "Voronoi Tessellation"
    @test TDGL_Polycrystal.state(f_2d_p) isa State
    @test TDGL_Polycrystal.system(f_2d_p) isa System

    TDGL_Polycrystal.jc_bcs!(f_2d_p, TDGL_Polycrystal.system(f_2d_p))
    @test all(Adapt.adapt(CPU(), data(f_2d_p.δda_rhs)) .<= 0)

    # 2D Non-periodic (JosephsonJunction) setup & BCs
    jj = TDGL_Polycrystal.JosephsonJunction(2.0, 6)
    f_2d_np, meta_2d_np, _, _, _, _ = simulation_setup(
        1, jj, tstep, GL, init_σ, norm_res, norm_mass, Ecrit, Jramp,
        wait_time, init_hold_time, 8, 16, false, alphaN, betaN, 1.0, 1.0,
        JcFinder, 2, tol, "CPU", 1, 1e-4, B_app, "LIN_INCREASE")
    TDGL_Polycrystal.jc_bcs!(f_2d_np, TDGL_Polycrystal.system(f_2d_np))
    @test sum(abs.(Adapt.adapt(CPU(), data(f_2d_np.δda_rhs)))) > 0

    # Unsupported periodicity branch
    bad_mesh = MulTDGL.RectMesh((8, 8), (false, false), (1.0, 1.0))
    bad_sys = System(bad_mesh, TDGL_Polycrystal.system(f_2d_np).p, TDGL_Polycrystal.system(f_2d_np).mat,
                     TDGL_Polycrystal.system(f_2d_np).c, TDGL_Polycrystal.system(f_2d_np).dc, TDGL_Polycrystal.system(f_2d_np).u)
    redirect_stdout(devnull) do
        TDGL_Polycrystal.jc_bcs!(f_2d_np, bad_sys)
    end

    # 3D Setup (zperiodic = false)
    vor3d = TDGL_Polycrystal.Voronoi(13.0, 1.0, 1)
    f_3d, meta_3d, weights_3d = redirect_stdout(devnull) do
        simulation_setup_3D(
            1, vor3d, tstep, GL, init_σ, norm_res, norm_mass, 1.5, Ecrit, Jramp,
            wait_time, init_hold_time, 16, 16, 16, true, false, alphaN, betaN, 1.0, 1.0,
            JcFinder, 2, tol, "CPU", 1, 1e-4, B_app, "EXP_INCREASE")
    end
    @test f_3d isa JcFinder
    @test meta_3d["Anisotropy"] == 1.5
    @test size(weights_3d) == (16, 16, 17)
    TDGL_Polycrystal.jc_bcs!(f_3d, TDGL_Polycrystal.system(f_3d))
    @test sum(abs.(Adapt.adapt(CPU(), data(f_3d.δda_rhs)))) > 0

    # 3D Setup (zperiodic = true) to cover 3D periodic set_fluxons!
    f_3d_p, _, weights_3d_p = redirect_stdout(devnull) do
        simulation_setup_3D(
            1, vor3d, tstep, GL, init_σ, norm_res, norm_mass, 1.0, Ecrit, Jramp,
            wait_time, init_hold_time, 16, 16, 16, true, true, alphaN, betaN, 1.0, 1.0,
            JcFinder, 2, tol, "CPU", 1, 1e-4, B_app, "EXP_INCREASE")
    end
    @test size(weights_3d_p) == (16, 16, 16)
    @test sum(Adapt.adapt(CPU(), data(TDGL_Polycrystal.system(f_3d_p).c))) > 0

    # new_finder & apply_B_field!
    f_new = TDGL_Polycrystal.new_finder(f_2d_p, JcFinder, Ecrit, wait_time, init_hold_time,
                                        Jramp, 1e-4, 0.2, tol, 2, "CPU", 2, "EXP_INCREASE")
    @test f_new.B_field ≈ 0.2
    apply_B_field!(f_new, 0.4)
    @test f_new.B_field ≈ 0.4
end
