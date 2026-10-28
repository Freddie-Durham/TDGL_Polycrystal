const corners = [[-0.5, 0.5], [0.5, 0.5], [0.5, -0.5], [-0.5, -0.5]]

struct JosephsonJunction 
    thick::Float64
    width::Int64
end

function append_metadata!(metadata::Dict, pattern::JosephsonJunction)
    metadata["Pattern"] = "Josephson Junction"
    metadata["Junction Thickness (ξ)"] = pattern.thick
    metadata["Junction Width (ξ)"] = pattern.width
end

"returns true if the point P is inside the rectangle defined by perpendicular edges AB and BC"
function inside_rectangle(AB, BC, AP, BP)
    return dot(AB, AP) >= 0 && dot(BC, BP) >= 0 && dot(AB, AB) >= dot(AB, AP) && dot(BC, BC) >= dot(BC, BP)
end

"perform binary search along a line to find intersection with box"
function search_intersect(cur_in, cur_out, AB, BC, box, iterations=10)
    new_pos = cur_in
    for _ in 1:iterations
        new_pos = (cur_in .+ cur_out)/2
        if inside_rectangle(AB, BC, new_pos .- box[1], new_pos .- box[2])
            cur_in = new_pos
        else
            cur_out = new_pos
        end
    end
    return new_pos
end

triangle_area(A, B, C) = 0.5 * rectangle_area(A, B, C)

"Area of rectangle ABCD where AB and BC are perpendicular edges"
rectangle_area(A, B, C) = abs((B[1] - A[1]) * (C[2] - A[2]) - (C[1] - A[1]) * (B[2] - A[2]))

"A and B are adjacent corners with right angles at A and B"
function quadrilateral_area(A, B, C, D)
    rect1 = rectangle_area(A, B, C)
    rect2 = rectangle_area(D, A, B)
    if rect1 > rect2
        return rect2 + (rect1 - rect2) / 2
    else
        return rect1 + (rect2 - rect1) / 2
    end
end

function area_inside_box(point,box)
    AB = box[2] .- box[1]
    BC = box[3] .- box[2]
    corner_bool = map(i->inside_rectangle(AB,BC,point .+i .- box[1],point .+i .- box[2]),corners)
    total = sum(corner_bool)

    if total == 4
        return 1.0

    elseif total == 0
        return 0.0

    #three corners are inside
    elseif total == 3
        ind = findfirst(c->!c,corner_bool)
        ind_left = mod1(ind - 1, 4)
        ind_right = mod1(ind + 1, 4)

        cur_out = point .+ corners[ind]
        cur_inR = point .+ corners[ind_right]
        cur_inL = point .+ corners[ind_left]

        insideR = search_intersect(cur_inR, cur_out, AB, BC, box)
        insideL = search_intersect(cur_inL, cur_out, AB, BC, box)

        return 1.0 - triangle_area(cur_out,insideR,insideL)

    #a single corner is inside
    elseif total == 1
        ind = findfirst(corner_bool)
        ind_left = mod1(ind - 1, 4)
        ind_right = mod1(ind + 1, 4)

        cur_in = point .+ corners[ind]
        cur_outR = point .+ corners[ind_right]
        cur_outL = point .+ corners[ind_left]

        insideR = search_intersect(cur_in, cur_outR, AB, BC, box)
        insideL = search_intersect(cur_in, cur_outL, AB, BC, box)

        return triangle_area(cur_in, insideR, insideL)

    #two corners are inside
    elseif total == 2
        insides = findall(corner_bool)
        outsides = findall(c->!c, corner_bool)

        if abs(insides[1] - outsides[1]) == 1
            cur_in1 = point .+ corners[insides[1]]
            cur_out1 = point .+ corners[outsides[1]]
            cur_in2 = point .+ corners[insides[2]]
            cur_out2 = point .+ corners[outsides[2]]
        else
            cur_in1 = point .+ corners[insides[1]]
            cur_out1 = point .+ corners[outsides[2]]
            cur_in2 = point .+ corners[insides[2]]
            cur_out2 = point .+ corners[outsides[1]]
        end

        inside1 = search_intersect(cur_in1, cur_out1, AB, BC, box)
        inside2 = search_intersect(cur_in2, cur_out2, AB, BC, box)

        return quadrilateral_area(cur_in1, cur_in2, inside2, inside1)
    end
end

function make_box(orig,dest,thick) 
    norm_dir = (dest - orig) ./ norm(dest - orig)
    perp_dir = [-norm_dir[2], norm_dir[1]]

    upper_start = orig .+ thick*perp_dir
    lower_start = orig .- thick*perp_dir

    upper_end = dest .+ thick*perp_dir
    lower_end = dest .- thick*perp_dir

    return [upper_start, lower_start, upper_end, lower_end]
end

function simple_line!(mesh, orig, dest, thick)
    dims = size(mesh)
    box = make_box(orig, dest, thick / 2)

    mins = floor.(Int32, reduce((a, b)->min.(a, b), box))
    maxes = ceil.(Int32, reduce((a, b)->max.(a, b), box))

    lims = [clamp(mins[i], 1, dims[i]):clamp(maxes[i], 1, dims[i]) for i in 1:2]
    boundingbox = CartesianIndices((lims[1], lims[2]))

    @inbounds @simd for I in boundingbox
        mesh[I] += area_inside_box(I.I, box)
    end
end

function apply_pattern(voronoi::Voronoi, pixels, veclen)
    package_dir = abspath(joinpath(@__DIR__, "..")) * "neper/2D_crystal/"
    return draw_lattice(veclen, package_dir, voronoi.seed, voronoi.grain_size * pixels, voronoi.grain_thick * pixels)
end

"Apply josephson junction pattern to a grid"
function apply_pattern(junc::JosephsonJunction, pixels, veclen)
    init_grid = zeros(Float64, veclen)

    junc_thick = Int(pixels * junc.thick)
    half_width = floor(Int, pixels * junc.width / 2)
    coat_width = Int(pixels * 5.0)

    # normal barrier
    init_grid[half_width+1:half_width+junc_thick, :] .= 1.0
    # normal lower coating
    init_grid[:, 1:coat_width] .= 1.0
    # normal upper coating
    init_grid[:, end-coat_width+1:end] .= 1.0
    return init_grid
end

"generate edges in the y direction, with extra set of edges for periodic boundary conditions"
function interp_edgesY(mesh, periodic_y)
    dims = size(mesh)
    new_mesh = zeros((dims[1], dims[2] - 1 + periodic_y))
    for x in 1:dims[1]
        for y in 2:dims[2]
            new_mesh[x, y-1] = (mesh[x, y-1] + mesh[x, y])/2
        end
        if periodic_y
            new_mesh[x, end] = (mesh[x, end] + mesh[x, 1])/2
        end
    end
    return new_mesh
end

"generate edges in the x direction, with extra set of edges for periodic boundary conditions"
function interp_edgesX(mesh, periodic_x)
    dims = size(mesh)
    new_mesh = zeros((dims[1] - 1 + periodic_x, dims[2]))
    for y in 1:dims[2]
        for x in 2:dims[1]
            new_mesh[x-1, y] = (mesh[x-1, y] + mesh[x, y])/2
        end
        if periodic_x
            new_mesh[end, y] = (mesh[end, y] + mesh[1, y])/2
        end
    end
    return new_mesh
end

struct Line{T <: AbstractVector}
    origin::T
    destination::T
end

line_inside_box(line::Line, dims) = (all(line.origin.>1) && all(line.origin.<dims)) || (all(line.destination.>1) && all(line.destination.<dims))

"process .tess file and draw anti-aliased lines in a periodic mesh. have to check for overdrawing of lines"
function create_2D_tessellation(filename, dims, GB_thick)
    lattice = TDGL_Polycrystal.Lattice2D(filename)
    mesh = zeros(dims)
    offsets = [[x, y] for x in [-1, 0, 1] * dims[1] for y in [-1, 0, 1] * dims[2]]

    line_dict = Dict{Tuple,Line}()
    for (_, edge) in lattice.edges
        orig = lattice.vertices[edge.p[1]].p .* dims
        dest = lattice.vertices[edge.p[2]].p .* dims

        for offset in offsets
            new_orig = orig .+ offset
            new_dest = dest .+ offset
            new_line = Line(new_orig, new_dest)

            if line_inside_box(new_line, dims)
                if new_orig[1] < new_dest[1]
                    approx_line = (round.(new_orig), round.(new_dest))
                else
                    approx_line = (round.(new_dest), round.(new_orig))
                end

                if !haskey(line_dict, approx_line)
                    line_dict[approx_line] = new_line
                    TDGL_Polycrystal.simple_line!(mesh, new_line.origin, new_line.destination, GB_thick)
                end
            end
        end
    end

    return min.(1.0, mesh) 
end

"takes in grain size and grain boundary thickness in pixels"
function draw_lattice(dims, filename, id, grain_size, GB_thick)
    grain_size /= dims[1]
    filename *= "D$(round(Int, grain_size * 1000))AR$(round(Int, 100 * dims[1]/dims[2]))xyID$id"
    if !isfile(filename*".tess")
        run(`neper -T -n from_morpho -morpho "graingrowth($grain_size)" -dim 2 -periodicity 1 -id $id -o $filename`)
    end
    return create_2D_tessellation(filename, dims, GB_thick)
end