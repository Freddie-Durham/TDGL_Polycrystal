@testset "Definitions & Mode Strings" begin
    d = Dict(:a => 42, :b => 99)
    @test d(:a) == 42
    @test string(JcInitHold()) == "Initial Hold"
    @test string(JcJHold()) == "Hold at fixed J"
    @test string(JcDone()) == "Simulation Complete"
end

@testset "Interpolation & Logging Utilities" begin
    @test TDGL_Polycrystal.linear_interp(10.0, 20.0, 0.5) ≈ 15.0
    @test TDGL_Polycrystal.linear_interp(-1.0, 1.0, 0.0) ≈ -1.0
    @test TDGL_Polycrystal.linear_interp(-1.0, 1.0, 1.0) ≈ 1.0

    p_empty = TDGL_Polycrystal.verbose_prefix("")
    @test p_empty isa String
    p_ctx = TDGL_Polycrystal.verbose_prefix("B=0.1")
    @test occursin("B=0.1", p_ctx)

    redirect_stdout(devnull) do
        TDGL_Polycrystal.verbose_update(["line 1", "line 2"]; context="test_ctx")
    end
end


@testset "HDF5 Metadata, Data & Checkpointing I/O" begin
    mktempdir() do tmpdir
        jj = TDGL_Polycrystal.JosephsonJunction(2.0, 6)
        f_io, meta, s_α, s_β, s_mx, s_my = simulation_setup(
            1, jj, 1.0, GL, init_σ, norm_res, norm_mass, Ecrit, Jramp,
            2, 1, 8, 16, false, alphaN, betaN, 1.0, 1.0,
            JcFinder, 2, tol, "CPU", 1, 1e-4, 0.25, "RAMP_B")

        path = tmpdir * "/"
        name = "test_sim"
        filepath = path * name * ".h5"
        chk_file = path * name * "_chk.h5"

        # save_metadata & save_simdata
        save_metadata(path, name, meta, s_α, s_β, s_mx, s_my)
        @test isfile(filepath)
        header = ["Current", "Electric Field", "Magnetic Field"]
        save_simdata(path, name, [[1e-4, 2e-4], [1e-6, 2e-6], [0.25, 0.25]], header, 1.23)

        # 3D key_metadata (weights) and appending in save_data
        HDF5.h5open(filepath, "r+") do fid
            g3d = HDF5.create_group(fid, "grid3d")
            TDGL_Polycrystal.key_metadata(g3d, ones(Float32, 4, 4, 4), Dict("Test" => 1))
            dg = fid["data"]
            TDGL_Polycrystal.save_data(header, [[3e-4], [3e-6], [0.25]], dg)
            @test TDGL_Polycrystal.saved_rows(header, dg) == 3
        end

        # write_checkpoint & save_state overwrite
        f_io.j = 5e-4
        f_io.B_field = 0.25
        f_io.mode = JcJHold()
        f_io.curholdsteps = 7
        f_io.esum = 1.5e-5
        f_io.E_field = 2.5e-6
        TDGL_Polycrystal.write_checkpoint(f_io, chk_file, 2)
        @test isfile(chk_file)
        HDF5.h5open(chk_file, "r+") do fid
            TDGL_Polycrystal.save_state(f_io, fid["save_state"], f_io.j, f_io.B_field)
        end

        # sync_data_with_checkpoint: truncate 3 rows back to 2 rows stored in checkpoint
        TDGL_Polycrystal.sync_data_with_checkpoint(filepath, chk_file, "data", header)
        HDF5.h5open(filepath, "r") do fid
            @test TDGL_Polycrystal.saved_rows(header, fid["data"]) == 2
        end
        # Warning branch when actual < expected
        TDGL_Polycrystal.write_checkpoint(f_io, chk_file, 5)
        @test_logs (:warn, r"fewer rows") TDGL_Polycrystal.sync_data_with_checkpoint(filepath, chk_file, "data", header)

        # Restore into a fresh finder
        f_restored = TDGL_Polycrystal.new_finder(f_io, JcFinder, Ecrit, 2, 1, Jramp, 1e-4, 0.25, tol, 2, "CPU", 2, "RAMP_B")
        HDF5.h5open(chk_file, "r") do fid
            TDGL_Polycrystal.load_state!(f_restored, fid["save_state"])
        end
        @test f_restored.j ≈ 5e-4
        @test f_restored.B_field ≈ 0.25
        @test f_restored.mode == JcJHold()
        @test f_restored.curholdsteps == 7
        @test f_restored.esum ≈ 1.5e-5
        @test f_restored.E_field ≈ 2.5e-6

        # Unknown mode warning and missing dataset error in load_state!
        HDF5.h5open(chk_file, "r+") do fid
            sg = fid["save_state"]
            HDF5.delete_attribute(sg, "Mode")
            HDF5.attributes(sg)["Mode"] = "UnknownMode"
        end
        HDF5.h5open(chk_file, "r") do fid
            @test_logs (:warn, r"Unknown mode") TDGL_Polycrystal.load_state!(f_restored, fid["save_state"])
        end
        HDF5.h5open(chk_file, "r+") do fid
            HDF5.delete_object(fid["save_state"], "ψ")
        end
        HDF5.h5open(chk_file, "r") do fid
            @test_throws ErrorException TDGL_Polycrystal.load_state!(f_restored, fid["save_state"])
        end
    end
end

