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
