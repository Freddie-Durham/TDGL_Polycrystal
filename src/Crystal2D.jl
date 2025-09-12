const corners = [[-0.5,0.5],[0.5,0.5],[0.5,-0.5],[-0.5,-0.5]]

#grain size and GB thickness in coherence lengths
struct Voronoi 
    grain_size::Float64
    grain_thick::Float64
    seed::Int64
end

function append_metadata!(metadata::Dict,pattern::Voronoi)
    metadata["Pattern"] = "Voronoi Tessellation"
    metadata["Grain size (ξ)"] = pattern.grain_size
    metadata["Grain Boundary Thickness (ξ)"] = pattern.grain_thick
end

#grain size and GB thickness in coherence lengths
struct TruncOct 
    grain_size::Float64
    grain_thick::Float64
    rep_grain::Int64
    angle::Float64
    factor::Int64
    shape::Any
end

function TruncOct(N::Integer,xmin::Integer,rep_grain::Integer,grain_thick::Float64,factor::Integer)
    grainangle,grain_diameter = periodic_crystal(N,fld(xmin,rep_grain))
    return TruncOct(grain_diameter,grain_thick,rep_grain,grainangle,factor,octagon!)
end

function append_metadata!(metadata::Dict,pattern::TruncOct)
    metadata["Pattern"] = "Truncated Octahedra"
    metadata["Grain size (ξ)"] = pattern.grain_size
    metadata["Grain Boundary Thickness (ξ)"] = pattern.grain_thick
    metadata["Lattice angle (rads)"] = pattern.angle
    metadata["Anti-aliasing factor"] = pattern.factor
    metadata["Multiple of grains"] = pattern.rep_grain
end

struct JosephsonJunction 
    junc_thick::Float64
end

"Returns size and angle of grain to enable periodic BCs"
function periodic_crystal(n,width)
    θ = atan(1/n)
    a = width/(n*cos(θ) + sin(θ))
    return θ, a
end

"returns true if the point P is inside the rectangle defined by perpendicular edges AB and BC"
function inside_rectangle(AB,BC,AP,BP)
    return dot(AB,AP) >= 0 && dot(BC,BP) >= 0 && dot(AB,AB) >= dot(AB,AP) && dot(BC,BC) >= dot(BC,BP)
end

"perform binary search along a line to find intersection with box"
function search_intersect(cur_in,cur_out,AB,BC,box,iterations=10)
    new_pos = cur_in
    for _ in 1:iterations
        new_pos = (cur_in .+ cur_out)/2
        if inside_rectangle(AB,BC,new_pos.-box[1],new_pos.-box[2])
            cur_in = new_pos
        else
            cur_out = new_pos
        end
    end
    return new_pos
end

triangle_area(A,B,C) = 0.5*rectangle_area(A,B,C)

"Area of rectangle ABCD where AB and BC are perpendicular edges"
rectangle_area(A,B,C) = abs((B[1]-A[1])*(C[2]-A[2])-(C[1]-A[1])*(B[2]-A[2]))

"A and B are adjacent corners with right angles at A and B"
function quadrilateral_area(A,B,C,D)
    rect1 = rectangle_area(A,B,C)
    rect2 = rectangle_area(D,A,B)
    if rect1 > rect2
        return rect2 + (rect1-rect2)/2
    else
        return rect1 + (rect2-rect1)/2
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
        ind_left = mod1(ind-1,4)
        ind_right = mod1(ind+1,4)

        cur_out = point .+ corners[ind]
        cur_inR = point .+ corners[ind_right]
        cur_inL = point .+ corners[ind_left]

        insideR = search_intersect(cur_inR,cur_out,AB,BC,box)
        insideL = search_intersect(cur_inL,cur_out,AB,BC,box)

        return 1.0 - triangle_area(cur_out,insideR,insideL)

    #a single corner is inside
    elseif total == 1
        ind = findfirst(corner_bool)
        ind_left = mod1(ind-1,4)
        ind_right = mod1(ind+1,4)

        cur_in = point .+ corners[ind]
        cur_outR = point .+ corners[ind_right]
        cur_outL = point .+ corners[ind_left]

        insideR = search_intersect(cur_in,cur_outR,AB,BC,box)
        insideL = search_intersect(cur_in,cur_outL,AB,BC,box)

        return triangle_area(cur_in,insideR,insideL)

    #two corners are inside
    elseif total == 2
        insides = findall(corner_bool)
        outsides = findall(c->!c,corner_bool)

        if abs(insides[1]-outsides[1]) == 1
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

        inside1 = search_intersect(cur_in1,cur_out1,AB,BC,box)
        inside2 = search_intersect(cur_in2,cur_out2,AB,BC,box)

        return quadrilateral_area(cur_in1,cur_in2,inside2,inside1)
    end
end

function make_box(orig,dest,thick) 
    norm_dir = (dest-orig)./norm(dest-orig)
    perp_dir = [-norm_dir[2],norm_dir[1]]

    upper_start = orig .+ thick*perp_dir
    lower_start = orig .- thick*perp_dir

    upper_end = dest .+ thick*perp_dir
    lower_end = dest .- thick*perp_dir

    return [upper_start,lower_start,upper_end,lower_end]
end

function simple_line!(mesh, orig, dest, thick)
    dims = size(mesh)

    box = make_box(orig,dest,thick/2)
    box_x = map(x->round(Int,x[1]),box)
    box_y = map(x->round(Int,x[2]),box)

    for xpxl in minimum(box_x):maximum(box_x)
        for ypxl in minimum(box_y):maximum(box_y)
            if is_safe([xpxl,ypxl],dims)
                mesh[xpxl,ypxl] += area_inside_box([xpxl,ypxl],box)
            end
        end
    end
end


"Calculate how close a point is to a 2D rectangular slab"
function isclose(x,y,dx,dy,thickness)
    if dy == 0
        xoff = 0
    else
        xoff = x - y*dx/dy
    end
    if dx==0
        yoff = 0
    else
        yoff = y - x*dy/dx
    end
    return (xoff^2 + yoff^2) < thickness^2
end

"Returns the endpoint of a line of certain length and angle from orig "
function get_dest(orig,length,angle)
    dx = Int(round(length*cos(angle)))
    dy = Int(round(length*sin(angle)))
    
    return [orig[1] + dx,orig[2] + dy]
end

"Draws a line in 2D in pixels"
function drawline!(grid,orig,dest,thickness=1,value=0)
    dx = dest[1] - orig[1]
    dy = dest[2] - orig[2]

    idir = dsign(dx)
    jdir = dsign(dy)

    #Deal with special cases of vertical and horizontal lines
    is_vert = dx==0 ? 1 : 0
    is_horiz = dy==0 ? 1 : 0
    domain_thickness = Int(ceil(thickness/2))

    diff = Int.(floor.(thickness/2(sqrt(dx^2+dy^2))*[dx,dy]))
    newdest = dest .- diff
    neworig = orig .+ diff

    imin = neworig[1] - is_vert * domain_thickness - idir*domain_thickness
    imax = newdest[1] + is_vert * domain_thickness + idir*domain_thickness
    jmin = neworig[2] - is_horiz * domain_thickness - jdir*domain_thickness
    jmax = newdest[2] + is_horiz * domain_thickness + jdir*domain_thickness

    for i in imin:idir:imax
        for j in jmin:jdir:jmax
            if isclose(i - neworig[1],j - neworig[2],dx,dy,thickness) && issafe(i,size(grid)[1]) && issafe(j,size(grid)[2])
                grid[i,j] = value
            end
        end
    end
end

"Draw a square using drawline!"
function square!(grid,origin,len,thickness,angle,value=0)
    dest = get_dest(origin,len,angle)
    drawline!(grid,origin,dest,thickness,value)

    for a in [π/2,-π,3*π/2]
        origin = dest
        dest = get_dest(origin,len,a+angle)
        drawline!(grid,origin,dest,thickness,value)
    end
end

"Draw an octagon using drawline!"
function octagon!(grid,origin,len,thickness,angle,value=0)
    dest = get_dest(origin,len/4,angle) #Starts from bottom right corner

    for a in [0,π/2,-π,3*π/2]
        origin = dest
        dest = get_dest(origin,len/2,angle+a)
        drawline!(grid,origin,dest,thickness,value)
        origin = dest
        dest = get_dest(origin,sqrt(2)*len/4,π/4+angle+a)
        drawline!(grid,origin,dest,thickness,value)
    end
end

"Draw a series of identical shapes that tesselate space"
function tesselate!(shape,grid,len,thickness,angle,xlen,ylen,value=0)
    safety_factor = 2 #min = sqrt(2)
    max_len = max(xlen,ylen)
    long_side = safety_factor*max_len
    next_i = [Int(max_len/2 - round(max_len*(1+safety_factor/sqrt(8))*cos(π/4+angle))),Int(max_len/2 - round(max_len*(1+safety_factor/sqrt(8))*sin(π/4+angle)))]

    for i in 1:Int(ceil(long_side/len))
        shape(grid,next_i,len,thickness,angle,value)
        next_j = get_dest(next_i,len,π/2 + angle)

        for j in 1:Int(ceil(long_side/len))
            shape(grid,next_j,len,thickness,angle,value)
            next_j = get_dest(next_j,len,π/2 + angle)
        end
        next_i = get_dest(next_i,len,angle)
    end 
end

"Applies non_periodic boundary conditions along top and bottom with lerp smoothing"
function non_periodic!(mesh,dims,ymin,ymax,startval=0,endval=1)
    for i in 1:dims[1]
        for j in ymin:ymax
            lrp = lerp(ymin,ymax,j)
            mesh[i,j] = startval*(1-lrp) + endval*lrp
            mesh[i,dims[2]+ymin-j] = startval*(1-lrp) + endval*lrp
        end
    end
end

"Apply tesselating pattern to a grid at high resolution then return anti-aliased scaled down version"
function apply_pattern(TrOct::TruncOct,init_val,new_val,pixels,veclen)
    init_grid = init_val*ones(Float64,veclen.*TrOct.factor)
    tesselate!(TrOct.shape,init_grid,TrOct.grain_size*pixels*TrOct.factor,TrOct.grain_thick*pixels*TrOct.factor,TrOct.angle,veclen[1]*TrOct.factor,veclen[2]*TrOct.factor,new_val)
    return lower_resolution(init_grid,TrOct.factor)
end

function apply_pattern(voronoi::Voronoi,init_val,new_val,pixels,veclen)
    package_dir = abspath(joinpath(@__DIR__,".."))
    return draw_lattice(veclen,package_dir*"neper/2D_crystal",voronoi.seed,voronoi.grain_size*pixels,voronoi.grain_thick*pixels,init_val,new_val)
end

"Apply josephson junction pattern to a grid"
function apply_pattern(JJ::JosephsonJunction,init_val,new_val,pixels,veclen)
    init_grid = init_val*ones(Float64,veclen)

    half_x = fld(veclen[1],2)
    junc_thick_xi = fld(JJ.junc_thick*pixels,2)
    init_grid[half_x-junc_thick_xi:half_x+junc_thick_xi,end - fld(veclen[2],3):end] .= new_val
    init_grid[half_x-junc_thick_xi:half_x+junc_thick_xi,1:fld(veclen[2],3)] .= new_val
    return init_grid
end

struct Vertex{K <: AbstractFloat} 
    p::Vector{K}
end
Vertex(;dim=2) = Vertex(zeros(Float64,dim))

struct Edge{K <: Integer}
    p::Vector{K}
end
Edge(;dim=2) = Edge(zeros(Int64,dim))

identifier(::Vertex) = "**vertex"
identifier(::Edge) = "**edge"

struct Lattice2D
    vertices::Dict{Int,Vertex}
    edges::Dict{Int,Edge}
end

function get_points(arr,type::Union{Vertex{K},Edge{K}}) where K
    ind = findfirst(s->contains(s,identifier(type)),arr)
    num_points = parse(Int64,arr[ind+1])

    T = typeof(type)
    dict = Dict{Int,T}()

    for i in ind+2:ind+1+num_points
        data = split(arr[i])
        p = T([parse(K,data[2]),parse(K,data[3])])
        dict[parse(Int,data[1])] = p
    end
    return dict
end

function Lattice2D(filename)
    strs = readlines(filename*".tess")
    Lattice2D(get_points(strs,Vertex()),get_points(strs,Edge()))
end

function draw!(mesh,orig,dest,GB_thick,val)
    o = map(x->round(Int,x),orig)
    d = map(x->round(Int,x),dest)

    if o!=d
        drawline!(mesh,o,d,GB_thick,val)
    end
end

"takes in grain size and grain boundary thickness in pixels"
function draw_lattice(dims,filename,id,grain_size,GB_thick,init_val=0.0,set_val=1.0)
    grain_size /= dims[1]
    filename *= "D$(round(Int,grain_size*1000))AR$(round(Int,100*dims[1]/dims[2]))xyID$id"
    if !isfile(filename*".tess")
        run(`neper -T -n from_morpho -morpho "graingrowth($grain_size)" -dim 2 -periodicity 1 -id $id -o $filename`)
    end

    lattice = Lattice2D(filename)
    mesh = zeros(dims)

    for (_,edge) in lattice.edges
        orig = lattice.vertices[edge.p[1]].p.*dims
        dest = lattice.vertices[edge.p[2]].p.*dims
        simple_line!(mesh, orig, dest, GB_thick)

        #apply periodicity condition
        if any(dest .> dims) || any(orig .> dims) || any(dest .< 1) || any(orig .< 1)
            orig1 = copy(orig)
            dest1 = copy(dest)
            if dest[1] > dims[1] || orig[1] > dims[1]
                orig1[1] -= dims[1]
                dest1[1] -= dims[1]
            end
            if dest[2] > dims[2] || orig[2] > dims[2]
                orig1[2] -= dims[2]
                dest1[2] -= dims[2]
            end
            if dest[1] < 1 || orig[1] < 1
                orig1[1] += dims[1]
                dest1[1] += dims[1]
            end
            if dest[2] < 1 || orig[2] < 1
                orig1[2] += dims[2]
                dest1[2] += dims[2]
            end
            simple_line!(mesh, orig1, dest1, GB_thick)
        end
    end

    return map((x)->linear_interp(init_val,set_val,min(1.0,x)),mesh) 
end