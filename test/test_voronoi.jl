@testset "Voronoi Metadata" begin
    meta = Dict{String, Any}()
    vor = TDGL_Polycrystal.Voronoi(13.0, 1.0, 1)
    TDGL_Polycrystal.append_metadata!(meta, vor)
    @test meta["Pattern"] == "Voronoi Tessellation"
    @test meta["Grain size (ξ)"] == 13.0
    @test meta["Grain Boundary Thickness (ξ)"] == 1.0
    @test meta["Voronoi Seed"] == 1
end

@testset "3D Triangle Geometry & Stencil" begin
    a = TDGL_Polycrystal.vec3(1.0, 1.0, 1.0)
    b = TDGL_Polycrystal.vec3(5.0, 1.0, 1.0)
    c = TDGL_Polycrystal.vec3(1.0, 5.0, 1.0)
    tdata = TDGL_Polycrystal.TriangleData(a, b, c, 2.0)

    @test TDGL_Polycrystal.inside_triangle(TDGL_Polycrystal.vec3(2.0, 2.0, 1.0), tdata) == true
    @test TDGL_Polycrystal.inside_triangle(TDGL_Polycrystal.vec3(2.0, 2.0, 5.0), tdata) == false
    @test TDGL_Polycrystal.inside_triangle(TDGL_Polycrystal.vec3(10.0, 10.0, 1.0), tdata) == false

    # Full pixel inside, boundary pixel (precise_area), and outside pixel
    @test TDGL_Polycrystal.triangle_weight(TDGL_Polycrystal.vec3(2.0, 2.0, 1.0), tdata) == 1.0f0
    @test 0.0f0 < TDGL_Polycrystal.triangle_weight(TDGL_Polycrystal.vec3(1.0, 1.0, 1.0), tdata) < 1.0f0
    @test TDGL_Polycrystal.triangle_weight(TDGL_Polycrystal.vec3(10.0, 10.0, 10.0), tdata) == 0.0f0
end

@testset "3D Edge Interpolation & Voronoi Mesh" begin
    m3d = ones(Float32, 4, 4, 4)
    mx, my, mz = TDGL_Polycrystal.interpolate_edges(m3d, (true, false, true))
    @test size(mx) == (4, 4, 4)
    @test size(my) == (4, 3, 4)
    @test size(mz) == (4, 4, 4)
    @test all(mx .≈ 1.0f0)

    # Uses existing cached neper/3D_crystal/D8125ID1.tess (13.0 * 1 / 16 = 0.8125)
    vor3d = TDGL_Polycrystal.Voronoi(13.0, 1.0, 1)
    w_vor3d = redirect_stdout(devnull) do
        TDGL_Polycrystal.apply_3D_pattern(vor3d, 1, (16, 16, 16))
    end
    @test size(w_vor3d) == (16, 16, 16)
    @test all(0.0f0 .<= w_vor3d .<= 1.0f0)
    @test any(w_vor3d .> 0.0f0)
end
