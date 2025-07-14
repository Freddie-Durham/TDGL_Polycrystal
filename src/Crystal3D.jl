"Returns whether a point is inside a triangular prism of given thickness"
function inside_triangle(p,a,b,c,thickness)::Bool
    vec_ab = b - a
    vec_ac = c - a
    tri_normal = normalize(cross(vec_ab,vec_ac))

    dist = dot(p,tri_normal) - dot(a,tri_normal) #signed perp dist to plane from point

    if dist^2 < (thickness^2)/4  #check if within thickness of plane
        p_proj = p - tri_normal*dist #find projection of point onto plane
        vec_ap = p_proj - a

        #compute barycentric coords
        abab = dot(vec_ab,vec_ab)
        acac = dot(vec_ac,vec_ac)
        abac = dot(vec_ab,vec_ac)
        apab = dot(vec_ap,vec_ab)
        apac = dot(vec_ap,vec_ac)

        inv_denom = 1/(abab * acac - abac * abac)
        u = (abab * apac - abac * apab)*inv_denom
        v = (acac * apab - abac * apac)*inv_denom

        if u>=0 && v>=0 && (v+u<=1) #check if within triangle inside plane
            return true
        else
            return false
        end
    else
        return false
    end
end

"Create a triangular prism of a given thickness"
function surface!(grid,a,b,c,thickness,value)
    xmax = Int(ceil(max(max(a[1],b[1]),c[1])))
    ymax = Int(ceil(max(max(a[2],b[2]),c[2])))
    zmax = Int(ceil(max(max(a[3],b[3]),c[3])))
    xmin = Int(floor(min(min(a[1],b[1]),c[1])))
    ymin = Int(floor(min(min(a[2],b[2]),c[2])))
    zmin = Int(floor(min(min(a[3],b[3]),c[3])))

    for i in xmin-thickness:xmax+thickness
        for j in ymin-thickness:ymax+thickness
            for k in zmin-thickness:zmax+thickness
                if issafe(i,size(grid)[1]) && issafe(j,size(grid)[2]) && issafe(k,size(grid)[3]) &&
                (grid[i,j,k]!=value) && inside_triangle([i,j,k],a,b,c,thickness)
                    grid[i,j,k] = value
                end
            end
        end
    end
end

"Check if in 3D grid"
function in_grid(val,upper,lower,limit)
    if all(val.+limit.>lower) && all(val.-limit.<upper)
        return true
    else
        return false
    end
end

"Create a square prism of a given thickness without rotating"
function face!(grid,start,basis,length,thickness,value)
    #
    #  3------4
    #  |      |
    #  |      |
    #  1------2
    # in x-y plane
    # basis defines vectors 1->2 and 1->3 respectively
    # length of side = length

    half = length/2

    p1 = start .- half*basis[1] .- half*basis[2]
    p2 = start .+ half*basis[1] .- half*basis[2]
    p3 = start .- half*basis[1] .+ half*basis[2]
    p4 = start .+ half*basis[1] .+ half*basis[2]

    surface!(grid,p1,p2,p4,thickness,value)
    surface!(grid,p1,p3,p4,thickness,value)
end

"create hexagon in a given basis space"
function hexagon!(grid,start,basis,length,thickness,value)
    #    5----6
    #   /      \
    #  /        \
    # 3-----7----4
    #  \        /
    #   \      /
    #    1----2
    # Across from 7 to 4 = length
    # Up from 7 to 5.5 = length * √3/2
    # Across from 5.5 to 6 = length / 2
    # start at 7, in x-y plane

    up = basis[2] * length * sqrt(3)/2
    right = basis[1] * length

    #go around anticlockwise
    points = [start + right,    #4                    
        start + right/2 - up,   #2             
        start - right/2 - up,   #1            
        start - right,          #3              
        start - right/2 + up,   #5              
        start + right/2 + up]   #6

    #create 6 small triangles from centre
    for i in 1:6
        j = i% 6 + 1
        surface!(grid,start,points[i],points[j],thickness,value)
    end
end

orthonorm_basis() = [[1,0,0],
                     [0,1,0],
                     [0,0,1]]

"Create a truncated octahedron in a given basis with size equal to distance between opposite square faces"
function trunc_oct!(grid,centre,basis,size,thickness,value)
    centre_to_square = size/2
    centre_to_hexagon = sqrt(3)*size/4

    #corresponds to both square and hexagon edge length
    edge_length = size/(2*sqrt(2))
    #Tilt from square to hexagon
    square_to_hexagon = 0.6959*π   
    
    face!(grid,centre+basis[3]*centre_to_square,basis,edge_length,thickness,value)
    face!(grid,centre-basis[3]*centre_to_square,basis,edge_length,thickness,value)
    
    newQ = attitude_quaternion(π/2,basis[1]+basis[2])
    newbasis = rotate_basis(basis,newQ)
    face!(grid,centre+newbasis[3]*centre_to_square,newbasis,edge_length,thickness,value)
    face!(grid,centre-newbasis[3]*centre_to_square,newbasis,edge_length,thickness,value)

    newQ = attitude_quaternion(π/2,basis[1]-basis[2])
    newbasis = rotate_basis(basis,newQ)
    face!(grid,centre+newbasis[3]*centre_to_square,newbasis,edge_length,thickness,value)
    face!(grid,centre-newbasis[3]*centre_to_square,newbasis,edge_length,thickness,value)
    
    newQ = attitude_quaternion(square_to_hexagon,basis[1])
    newbasis = rotate_basis(basis,newQ)
    hexagon!(grid,centre+newbasis[3]*centre_to_hexagon,newbasis,edge_length,thickness,value)
    hexagon!(grid,centre-newbasis[3]*centre_to_hexagon,newbasis,edge_length,thickness,value)

    newQ = attitude_quaternion(-square_to_hexagon,basis[1])
    newbasis = rotate_basis(basis,newQ)
    hexagon!(grid,centre+newbasis[3]*centre_to_hexagon,newbasis,edge_length,thickness,value)
    hexagon!(grid,centre-newbasis[3]*centre_to_hexagon,newbasis,edge_length,thickness,value)
    
    newQ = attitude_quaternion(square_to_hexagon,basis[2])
    newbasis = rotate_basis(basis,newQ)
    hexagon!(grid,centre+newbasis[3]*centre_to_hexagon,[newbasis[2],newbasis[1]],edge_length,thickness,value)
    hexagon!(grid,centre-newbasis[3]*centre_to_hexagon,[newbasis[2],newbasis[1]],edge_length,thickness,value)

    newQ = attitude_quaternion(-square_to_hexagon,basis[2])
    newbasis = rotate_basis(basis,newQ)
    hexagon!(grid,centre+newbasis[3]*centre_to_hexagon,[newbasis[2],newbasis[1]],edge_length,thickness,value)
    hexagon!(grid,centre-newbasis[3]*centre_to_hexagon,[newbasis[2],newbasis[1]],edge_length,thickness,value)
end

"Create a cube in a given basis space"
function cube!(grid,centre,basis,edge_len,thickness,value)
    """
    Start from centre, must place faces in 6 directions:
    +x - face in y,z
    -x - face in y,z
    +y - face in x,z
    -y - face in x,z
    +z - face in x,y
    -z - face in x,y
    These x,y,z are defined by orthonormal basis
    """
    x = basis[1]
    y = basis[2]
    z = basis[3]
    
    centre_to_face = edge_len/2

    for (dir,face_basis) in zip([x,y,z],[[y,z],[x,z],[x,y]])
        for sgn in [+1,-1]
            c = centre + sgn*centre_to_face * dir
            face!(grid,c,face_basis,edge_len,thickness,value)
        end
    end
end

#See bcc lattice points defining trunc octahedra
function visualise_BCC!(grid,pos,dims,value)
    if in_grid(pos,dims,[0,0,0],0)
        pos = Int.(round.(pos))
        grid[pos[1],pos[2],pos[3]] = value
    end
end

"Rotate and tesselate 3D object throughout domain"
function tesselateOct!(grid,dims,q::Quaternion,size,thickness=1,value=2)
    start_pos = dims./2

    #create linearly independant basis to draw a trunc-oct
    basis = rotate_basis(orthonorm_basis(),q)
    #need to rotate 45 degrees to fit inside BCC structure
    rotate_TO = attitude_quaternion(π/4,basis[3])
    basis = rotate_basis(basis,rotate_TO)

    max_cells = Int.(2*ceil.(dims./size))

    #primitive BCC lattice vectors to tesselate space
    xbcc = (size/2)*[-1,1,1]
    ybcc = (size/2)*[1,-1,1]
    zbcc = (size/2)*[1,1,-1]

    #rotate BCC lattice vectors 
    xbcc_rot = rotate(xbcc,q)
    ybcc_rot = rotate(ybcc,q)
    zbcc_rot = rotate(zbcc,q)

    #i denotes the current number of the cell in a given direction
    for i in -max_cells[1]:max_cells[1]
        for j in -max_cells[2]:max_cells[2]
            for k in -max_cells[3]:max_cells[3]
                pos = start_pos + 1*(i*xbcc_rot + j*ybcc_rot + k*zbcc_rot)
                #visualise_BCC!(grid,pos,dims,value)
                if in_grid(pos,dims,[0,0,0],size)
                    trunc_oct!(grid,pos,basis,size,thickness,value)
                end
            end
        end
    end
end

"Rotate and tesselate 3D object throughout domain"
function tesselateCube!(grid,dims,q::Quaternion,edge_len,thickness=1,value=2)
    start_pos = [0.0,0.0,0.0]

    #create linearly independant basis
    basis = [rotate([1,0,0],q),
             rotate([0,1,0],q),
             rotate([0,0,1],q)]

    max_cells = Int.(2*ceil.(dims./edge_len))
 
    #i denotes the current number of the cell in a given direction
    for i in -max_cells[1]:max_cells[1]
        for j in -max_cells[2]:max_cells[2]
            for k in -max_cells[3]:max_cells[3]
                pos = start_pos + edge_len*(i*basis[1] + j*basis[2] + k*basis[3])
                if in_grid(pos,dims,[0,0,0],edge_len)
                    cube!(grid,pos,basis,edge_len,thickness,value)
                end
            end
        end
    end
end

