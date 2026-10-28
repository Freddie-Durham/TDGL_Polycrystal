@testset "Command Line Parsing" begin
    old_args = copy(ARGS)
    empty!(ARGS)
    parsed = try
        parse_CL()
    finally
        append!(ARGS, old_args)
    end
    @test parsed[:dims] == 2
    @test parsed[:ramp_mode] == "EXP_INCREASE"
    @test parsed[:backend] == "CPU"
end