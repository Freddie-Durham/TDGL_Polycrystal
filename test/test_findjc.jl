@testset "Electric Field Calculation & Stepping" begin
    jj = TDGL_Polycrystal.JosephsonJunction(2.0, 6)
    f_2d, _, _, _, _, _ = simulation_setup(
        1, jj, 1.0, GL, init_σ, norm_res, norm_mass, 1e-3, 0.5,
        2, 1, 8, 16, false, alphaN, betaN, 1.0, 1.0,
        JcFinder, 2, tol, "CPU", 1, 1e-4, 0.1, "EXP_INCREASE")

    ef = TDGL_Polycrystal.E_field_avg(f_2d.solver, TDGL_Polycrystal.state(f_2d).a, 4)
    @test ef isa Float64

    TDGL_Polycrystal.step!(f_2d)
    @test f_2d.curholdsteps == 0
    @test f_2d.mode == JcJHold()
end

@testset "Ramp Modes & find_jc Execution" begin
    @test TDGL_Polycrystal.decode_method("EXP_INCREASE", 0.5) isa TDGL_Polycrystal.Exp_Increase
    @test TDGL_Polycrystal.decode_method("EXP_DECREASE", 0.5) isa TDGL_Polycrystal.Exp_Decrease
    @test TDGL_Polycrystal.decode_method("LIN_INCREASE", 0.5) isa TDGL_Polycrystal.Linear_Increase
    @test TDGL_Polycrystal.decode_method("RAMP_B", 0.5) isa TDGL_Polycrystal.Ramp_B
    @test_throws ErrorException TDGL_Polycrystal.decode_method("INVALID", 0.5)

    jj = TDGL_Polycrystal.JosephsonJunction(2.0, 6)
    f_base, _, _, _, _, _ = simulation_setup(
        1, jj, 1.0, GL, init_σ, norm_res, norm_mass, 1e-3, 0.5,
        2, 1, 8, 16, false, alphaN, betaN, 1.0, 1.0,
        JcFinder, 2, tol, "CPU", 1, 1e-4, 0.1, "EXP_INCREASE")

    # 1. Exp_Increase
    f_exp = TDGL_Polycrystal.new_finder(f_base, JcFinder, 10.0, 2, 1, 0.5, 1e-3, 0.1, tol, 2, "CPU", 1, "EXP_INCREASE")
    f_exp.mode = JcJHold()
    f_exp.E_field = 0.0
    TDGL_Polycrystal.next_step!(f_exp, f_exp.ramp_method, TDGL_Polycrystal.system(f_exp).p)
    @test f_exp.j ≈ 1e-3 * 1.5
    f_exp.E_field = 20.0
    f_exp.curholdsteps = 2
    TDGL_Polycrystal.next_step!(f_exp, f_exp.ramp_method, TDGL_Polycrystal.system(f_exp).p)
    @test f_exp.mode == JcDone()

    # 2. Linear_Increase
    f_lin = TDGL_Polycrystal.new_finder(f_base, JcFinder, 10.0, 2, 1, 0.005, 1e-3, 0.1, tol, 2, "CPU", 1, "LIN_INCREASE")
    f_lin.mode = JcJHold()
    f_lin.E_field = 0.0
    TDGL_Polycrystal.next_step!(f_lin, f_lin.ramp_method, TDGL_Polycrystal.system(f_lin).p)
    @test f_lin.j ≈ 1e-3 + 0.005
    f_lin.E_field = 20.0
    f_lin.curholdsteps = 2
    TDGL_Polycrystal.next_step!(f_lin, f_lin.ramp_method, TDGL_Polycrystal.system(f_lin).p)
    @test f_lin.mode == JcDone()

    # 3. Exp_Decrease
    f_dec = TDGL_Polycrystal.new_finder(f_base, JcFinder, 1.0, 4, 1, 0.2, 1e-2, 0.1, tol, 2, "CPU", 1, "EXP_DECREASE")
    f_dec.mode = JcJHold()
    f_dec.curholdsteps = 4
    f_dec.esum = 10.0
    f_dec.E_field = 5.0
    TDGL_Polycrystal.next_step!(f_dec, f_dec.ramp_method, TDGL_Polycrystal.system(f_dec).p)
    @test f_dec.j ≈ 1e-2 * 0.8
    @test f_dec.mode == JcJHold()
    f_dec.curholdsteps = 4
    f_dec.esum = 0.1
    f_dec.E_field = 0.1
    TDGL_Polycrystal.next_step!(f_dec, f_dec.ramp_method, TDGL_Polycrystal.system(f_dec).p)
    @test f_dec.mode == JcDone()

    # 4. Ramp_B
    f_rb = TDGL_Polycrystal.new_finder(f_base, JcFinder, 1.0, 2, 1, 0.5, 1e-3, 1.5e-3, tol, 2, "CPU", 1, "RAMP_B")
    @test f_rb.B_field == 0.0
    f_rb.mode = JcJHold()
    TDGL_Polycrystal.next_step!(f_rb, f_rb.ramp_method, TDGL_Polycrystal.system(f_rb).p)
    @test f_rb.B_field ≈ 1e-3
    @test f_rb.mode == JcJHold()
    TDGL_Polycrystal.next_step!(f_rb, f_rb.ramp_method, TDGL_Polycrystal.system(f_rb).p)
    @test f_rb.B_field ≈ 2e-3
    @test f_rb.mode == JcInitHold()
    f_rb.mode = JcJHold()
    f_rb.E_field = 0.0
    TDGL_Polycrystal.next_step!(f_rb, f_rb.ramp_method, TDGL_Polycrystal.system(f_rb).p)
    @test f_rb.j ≈ 1e-3 * 1.5
    f_rb.E_field = 5.0
    f_rb.curholdsteps = 2
    TDGL_Polycrystal.next_step!(f_rb, f_rb.ramp_method, TDGL_Polycrystal.system(f_rb).p)
    @test f_rb.mode == JcDone()

    # 5. Full find_jc run with checkpointing & verbose output
    mktempdir() do tmpdir
        f_run = TDGL_Polycrystal.new_finder(f_base, JcFinder, -1.0, 1, 1, 0.1, 1e-3, 0.1, tol, 2, "CPU", 1, "EXP_INCREASE")
        data_file = joinpath(tmpdir, "run_data.h5")
        chk_file = joinpath(tmpdir, "run_chk.h5")
        header = ["Current", "Electric Field", "Magnetic Field"]
        HDF5.h5open(data_file, "w") do fid
            g = HDF5.create_group(fid, "data")
            HDF5.create_group(g, "0.1b data")
        end
        sim_data, timetaken = redirect_stdout(devnull) do
            find_jc(f_run, true; filepath=data_file, checkpoint_file=chk_file,
                    save_states=true, save_frequency=1, B=0.1, header=header,
                    verbose_context="test", report_period=2)
        end
        @test f_run.mode == JcDone()
        @test isfile(chk_file)
        @test timetaken >= 0.0
    end
end
