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
    metadata["Voronoi Seed"] = pattern.seed
end

struct Vertex{N,K<:AbstractFloat} 
    p::SVector{N,K}
end
Vertex{N}(::Type{K}=Float64) where {N, K<:AbstractFloat} = Vertex(zeros(SVector{N,K}))
function (::Type{Vertex{N,K}})(data,_) where {N, K<:AbstractFloat} 
    return Vertex(SVector{N,K}(parse.(K, data[2:N+1])))
end

struct Edge{K<:Integer}
    p::Vector{K}
end
Edge() = Edge(zeros(Int64,2))
function (::Type{Edge{K}})(data,_) where K<:Integer
    return Edge(parse.(K,data[2:3]))
end

struct Face{K<:Integer}
    p::Vector{K}
end
Face() = Face(zeros(Int64,3))
function (::Type{Face{K}})(_,data) where K<:Integer
    return Face(parse.(K,data[2:end]))
end

identifier(::Vertex) = "**vertex"
identifier(::Edge) = "**edge"
identifier(::Face) = "**face"

struct Lattice3D{T<:AbstractFloat,K<:Integer}
    vertices::Dict{Int,Vertex{3,T}}
    edges::Dict{Int,Edge{K}}
    faces::Dict{Int,Face{K}}
end

struct Lattice2D
    vertices::Dict{Int,Vertex}
    edges::Dict{Int,Edge}
end

data_range(::Face,ind,num_points) = ind+2:4:ind+1+4*num_points
data_range(::Union{Vertex,Edge},ind,num_points) = ind+2:ind+1+num_points

"retrieve elements of degree 0, 1 or 2 from an array of lines from a .tess file"
function get_points(arr,type::Union{Vertex{K},Edge{K},Face{K}}) where {K}
    ind = findfirst(s->contains(s,identifier(type)),arr)
    num_points = parse(Int64,arr[ind+1])

    T = typeof(type)
    dict = Dict{Int,T}()

    for i in data_range(type,ind,num_points)
        data = split(arr[i])
        dict[parse(Int,data[1])] = T(data,split(arr[i+1]))
    end
    return dict
end

"read a .tess file and construct a 3D voronoi tessellation"
function Lattice3D(filename)
    strs = readlines(filename*".tess")
    Lattice3D(get_points(strs,Vertex{3}()),get_points(strs,Edge()),get_points(strs,Face()))
end

function Lattice2D(filename)
    strs = readlines(filename*".tess")
    Lattice2D(get_points(strs,Vertex{2}()),get_points(strs,Edge()))
end

struct TriangleData
    a::vec3
    normal::vec3
    vec_ab::vec3
    vec_ac::vec3
    abab::Float64
    acac::Float64
    abac::Float64
    inv_denom::Float64
    half_thick2::Float64
    plane_offset::Float64
end

"precompute data for triangle calculation"
function TriangleData(a,b,c,thickness)
    vec_ab = b - a
    vec_ac = c - a
    tri_normal = normalize(cross(vec_ab,vec_ac))
    abab = dot(vec_ab,vec_ab)
    acac = dot(vec_ac,vec_ac)
    abac = dot(vec_ab,vec_ac)
    inv_denom = 1/(abab * acac - abac * abac)
    plane_offset = dot(a,tri_normal)
    half_thick2 = (thickness^2)/4
    return TriangleData(a,tri_normal,vec_ab,vec_ac,
    abab,acac,abac,inv_denom,half_thick2,plane_offset)
end

"Returns true if a point p is inside a triangular prism (a,b,c) of given thickness"
function inside_triangle(p,T::TriangleData)::Bool
    dist = dot(p,T.normal) - T.plane_offset #signed perp dist to plane from point

    if dist^2 > T.half_thick2  #check if within thickness of plane
        return false
    end

    p_proj = p - T.normal*dist #find projection of point onto plane
    vec_ap = p_proj - T.a

    #compute barycentric coords
    apab = dot(vec_ap,T.vec_ab)
    apac = dot(vec_ap,T.vec_ac)

    u = (T.abab * apac - T.abac * apab)*T.inv_denom
    v = (T.acac * apab - T.abac * apac)*T.inv_denom

    if (u>=0) & (v>=0) & (v+u<=1) #check if within triangle inside plane
        return true
    else
        return false
    end
end

#All 8 corners of the cube
OFFSETS::Vector{vec3} = [vec3([0.5,0.5,0.5]),vec3([-0.5,-0.5,-0.5]),
vec3([-0.5,0.5,0.5]),vec3([0.5,-0.5,0.5]),vec3([0.5,0.5,-0.5]),
vec3([-0.5,0.5,-0.5]),vec3([-0.5,-0.5,0.5]),vec3([0.5,-0.5,-0.5])]

#NxNxN points distributed evenly throughout the cube 
const N = 5
const VOLUME_FRAC = 1/(N^3)
PRECISE::Vector{vec3} = [vec3(([i,j,k]/(N-1)).-0.5) for i in 0:N-1 for j in 0:N-1 for k in 0:N-1]

"Use a more expensive stencil to estimate volume of pixel inside triangle"
function precise_area(point,tri_data::TriangleData)
    val = 0.0
    for vec in PRECISE
        if inside_triangle(point.+vec,tri_data)
            val += VOLUME_FRAC
        end
    end
    return val
end

"return a value between 0 and 1 based on volume of pixel inside triangle"
function triangle_weight(point,tri_data::TriangleData)::Float64
    val = 0.0
    @simd for off in OFFSETS
        if inside_triangle(point.+off,tri_data)
            val += 0.125
        end
    end
    if val > 0.0 && val < 1.0
        return precise_area(point,tri_data)
    else
        return val
    end
end

"draw an anti-aliased triangular prism of a given thickness on a 3D mesh"
function draw_slab!(mesh,a,b,c,thick = 1)
    dims = size(mesh)
    maxes = ceil.(Int32,max.(max.(a,b),c) .+ thick)
    mins = floor.(Int32,min.(min.(a,b),c) .- thick)

    lims = [clamp(mins[i],1,dims[i]):clamp(maxes[i],1,dims[i]) for i in 1:3]
    boundingbox = CartesianIndices((lims[1],lims[2],lims[3]))
    tri_data = TriangleData(a,b,c,thick)
    
    @inbounds for I in boundingbox
        mesh[I] += triangle_weight(vec3(Tuple(I)),tri_data)
    end
end

"split a convex polygon into triangles"
function fan_triangle(verts)
    start = verts[1]
    triangles = Tuple[]
    for i in 2:length(verts)-1
        a = verts[i]
        b = verts[i+1]
        push!(triangles,(a,b,start))
    end
    return triangles
end

"take a face from a voronoi tessellation and return a list of vertices of triangles to construct the face"
function triangulate(tess::Lattice3D{T,K},face::Face) where {T,K}
    edges::Vector{Vector{K}} = map(face.p) do e
        edge = tess.edges(abs(e)).p
        return sign(e)>0 ? edge : [edge[2],edge[1]]
    end
    verts::Vector{Vector{T}} = map(e->tess.vertices[e[1]].p,edges)
    return fan_triangle(verts)
end

"loop over all 27 unit cells and attempt to draw all triangles. 
this is fine since when drawing triangles we calculate our bounding box to be strictly inside the mesh."
function draw_periodic_triangle!(mesh,a,b,c,thick)
    dims = size(mesh)
    shifts = [-1, 0, 1]

    for dx in shifts, dy in shifts, dz in shifts
        shift = (dx*dims[1], dy*dims[2], dz*dims[3])
        A,B,C = map(v->v.+shift,(a,b,c))
        
        draw_slab!(mesh,A,B,C,thick)
    end
end

tri_area(a,b,c) = 0.5 * norm(cross(b - a, c - a))

function create_3D_mesh(tess::Lattice3D,dims,thick,verbose=false)
    mesh = zeros(dims)

    start_t = time()
    actual_volume = 0
    for (_,face) in tess.faces
        triangles = triangulate(tess,face)
        for tri in triangles
            a,b,c = map(t->t.*dims,tri)

            draw_periodic_triangle!(mesh,a,b,c,thick)
            actual_volume += tri_area(a,b,c)*thick
        end
    end
    calc_volume = sum(mesh)
    δt = time() - start_t

    println("Created mesh with fractional GB volume = $(round(Int,actual_volume/prod(dims))) in $(round(δt,sigdigits=4))s")
    println("Percent Volume Error = $(round(100*(actual_volume-calc_volume)/actual_volume,sigdigits=4))%")
    return map(m->min(m,1.0),mesh)
end

function interpolate_edges(mesh,periodic)
    dims = size(mesh)
    meshx = zeros((dims[1]-1+periodic[1],dims[2],dims[3]))
    meshy = zeros((dims[1],dims[2]-1+periodic[2],dims[3]))
    meshz = zeros((dims[1],dims[2],dims[3]-1+periodic[3]))

    for I in CartesianIndices(meshx)
        meshx[I] = 0.5 * (mesh[I] + mesh[mod1(I[1]+1,dims[1]),I[2],I[3]])
    end

    for I in CartesianIndices(meshy)
        meshy[I] = 0.5 * (mesh[I] + mesh[I[1],mod1(I[2]+1,dims[2]),I[3]])
    end

    for I in CartesianIndices(meshz)
        meshz[I] = 0.5 * (mesh[I] + mesh[I[1],I[2],mod1(I[3]+1,dims[3])])
    end

    return meshx, meshy, meshz
end

"Call neper to create a 3D voronoi tessellation if it does not already exist and return filename"
function create_3D_tessellation(filepath,grain_size,seed)
    filename = filepath * "D$(round(Int,grain_size*10000))ID$(seed)"
    if !isfile(filename*".tess")
        run(`neper -T -n from_morpho -morpho "graingrowth($grain_size)" -dim 3 -periodicity 1 -id $seed -o $filename`)
    end
    return filename
end

"generate a 3D voronoi tessellation pattern using neper and apply it to a mesh"
function apply_3D_pattern(voronoi::Voronoi,pixels,dims)
    filepath = abspath(joinpath(@__DIR__,".."))*"neper/3D_crystal/"
    grain_size = voronoi.grain_size * pixels / dims[1]
    lattice  = Lattice3D(create_3D_tessellation(filepath,grain_size,voronoi.seed))

    println("Generating 3D Voronoi pattern with grain size $(voronoi.grain_size)ξ and thickness $(voronoi.grain_thick)ξ. Dims = $dims")
    return create_3D_mesh(lattice, dims, voronoi.grain_thick*pixels, true)
end