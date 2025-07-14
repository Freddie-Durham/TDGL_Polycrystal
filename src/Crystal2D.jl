"Returns size and angle of grain to enable periodic BCs"
function periodic_crystal(n,width)
    θ = atan(1/n)
    a = width/(n*cos(θ) + sin(θ))
    return θ, a
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

"Draw a series of identical shapes that tesselate space"
function tesselateNEW!(shape,grid,dims,len,thickness,angle,value=0)
    safety_factor = 2 #min = sqrt(2)

    startpos = Int.(dims/2 .- len/2)
    dirx = [cos(angle),-sin(angle)]
    diry = [sin(angle),cos(angle)]

    max_len = max(dims[1],dims[2])
    long_side = safety_factor*max_len
    intlen = Int(ceil(long_side/len))

    for i in -intlen:intlen
        for j in -intlen:intlen
            pos = Int.(startpos + i*dirx + j*diry)
            if in_grid(pos,dims,[0,0],len)
                shape(grid,pos,len,thickness,angle,value)
            end
        end
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