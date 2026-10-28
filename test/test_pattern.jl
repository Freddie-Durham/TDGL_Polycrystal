@testset "Box & Area Calculations" begin
    A = [0.0, 0.0]
    B = [2.0, 0.0]
    C = [2.0, 3.0]
    D = [0.0, 3.0]
    @test TDGL_Polycrystal.rectangle_area(A, B, C) ≈ 6.0
    @test TDGL_Polycrystal.triangle_area(A, B, C) ≈ 3.0
    @test TDGL_Polycrystal.quadrilateral_area(A, B, C, D) ≈ 6.0
    @test TDGL_Polycrystal.quadrilateral_area(A, B, [2.0, 4.0], [0.0, 2.0]) ≈ 6.0
    @test TDGL_Polycrystal.quadrilateral_area(A, B, [2.0, 2.0], [0.0, 4.0]) ≈ 6.0

    # make_box and area_inside_box branches (0, 1, 2, 3, 4 corners inside)
    box = TDGL_Polycrystal.make_box([0.0, 0.0], [10.0, 0.0], 2.0)
    @test length(box) == 4
    # 4 corners inside
    @test TDGL_Polycrystal.area_inside_box((5.0, 0.0), box) ≈ 1.0
    # 0 corners inside
    @test TDGL_Polycrystal.area_inside_box((5.0, 10.0), box) ≈ 0.0
    # 2 corners inside (straddling top edge y=2.0)
    a2 = TDGL_Polycrystal.area_inside_box((5.0, 2.0), box)
    @test 0.0 < a2 < 1.0
    # 2 corners inside (straddling left edge x=0.0, hitting the else branch of abs(insides[1]-outsides[1])==1)
    a2_left = TDGL_Polycrystal.area_inside_box((-0.2, 0.0), box)
    @test 0.0 < a2_left < 1.0
    # 1 corner inside (near corner (10.0, 2.0))
    a1 = TDGL_Polycrystal.area_inside_box((10.2, 2.2), box)
    @test 0.0 < a1 < 0.5
    # 3 corners inside (angled box cutting one corner of pixel)
    mesh = zeros(12, 12)
    TDGL_Polycrystal.simple_line!(mesh, [1.0, 1.0], [10.0, 8.0], 3.0)
    @test sum(mesh) > 0.0
end

@testset "2D Edge Interpolation" begin
    m2d = [1.0 3.0; 5.0 7.0]
    ex_p = TDGL_Polycrystal.interp_edgesX(m2d, true)
    ex_np = TDGL_Polycrystal.interp_edgesX(m2d, false)
    @test size(ex_p) == (2, 2)
    @test size(ex_np) == (1, 2)
    @test ex_np[1, 1] ≈ 3.0

    ey_p = TDGL_Polycrystal.interp_edgesY(m2d, true)
    ey_np = TDGL_Polycrystal.interp_edgesY(m2d, false)
    @test size(ey_p) == (2, 2)
    @test size(ey_np) == (2, 1)
    @test ey_np[1, 1] ≈ 2.0
end

@testset "2D Patterns: JosephsonJunction & Voronoi" begin
    meta = Dict{String, Any}()

    # JosephsonJunction
    jj = TDGL_Polycrystal.JosephsonJunction(2.0, 10)
    TDGL_Polycrystal.append_metadata!(meta, jj)
    @test meta["Pattern"] == "Josephson Junction"
    w_jj = apply_pattern(jj, 1, (12, 20))
    @test size(w_jj) == (12, 20)
    @test all(w_jj[:, 1:5] .== 1.0)
    @test all(w_jj[6:7, :] .== 1.0)
    @test w_jj[2, 10] == 0.0

    # Voronoi 2D (using existing D500AR100xyID1.tess: grain_size * pixels / dims[1] = 0.5)
    vor2d = TDGL_Polycrystal.Voronoi(8.0, 1.0, 1)
    w_vor2d = apply_pattern(vor2d, 1, (16, 16))
    @test size(w_vor2d) == (16, 16)
    @test all(0.0 .<= w_vor2d .<= 1.0)
    @test any(w_vor2d .> 0.0)
end
