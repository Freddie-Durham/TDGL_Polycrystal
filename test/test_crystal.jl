@testset "Periodic Cyrstal" begin
    θ, a = TDGL_Polycrystal.periodic_crystal(2,64)
    @test θ ≈ atan(0.5)
    @test a ≈ 28.621670111997307
end

@testset "Get Destination" begin
    dest = TDGL_Polycrystal.get_dest([0, 0], 5, π/4)
    @test dest == [4, 4]
end

@testset "Is Close" begin
    @test TDGL_Polycrystal.isclose(0.5, 0.5, 1, 1, 0.5) == true
    @test TDGL_Polycrystal.isclose(1, 2, 1, 0, 0.5) == false
end

@testset "Draw Line" begin
    grid = zeros(Int, 10, 10)
    TDGL_Polycrystal.drawline!(grid, [1, 1], [1, 5], 1, 1)
    @test sum(grid[1:5, 1]) == 1
    @test sum(grid[6:end, 1]) == 0
end