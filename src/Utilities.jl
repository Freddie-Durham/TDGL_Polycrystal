

"Used to rotate in 3D. Real number is angle in radians, imaginary 3Vec is axis of rotation"
struct Quaternion
    re::Real
    im::Vector{Real}
end

"Normalises quaternion and converts real value to rotation angle"
function attitude_quaternion(α,v::Vector{T}) where {T<:Real} #returns a normalised quaternion
    v = normalize(v)
    return Quaternion(cos(α/2),v*sin(α/2))
end

"Returns the complex conjugate of the quaternion"
conjugate(q::Quaternion) = Quaternion(q.re,-q.im) #Assumes normalised quaternion
"Returns absolute magnitude of quaternion"
magnitude_sq(q::Quaternion) = q.re^2 + q.im[1]^2 + q.im[2]^2 + q.im[3]^2

Base.:*(q::Quaternion,p::Quaternion) = Quaternion(q.re*p.re - dot(q.im,p.im),q.re*p.im + p.re*q.im + cross(q.im,p.im))

Base.show(q::Quaternion) = print("Real:",q.re," Imaginary:",q.im)

"Multiplies a normalised quaternion by a vector to return a vector"
function rotate(v,q::Quaternion) 
    v_quat = Quaternion(0,v)
    return (q*v_quat*conjugate(q)).im
end

"Rotate point p about r_point around r_vec with angle α"
function rotate_about_point(p,r_point,q::Quaternion) 
    return rotate(p-r_point,q) + r_point
end

"rotates a set of basis vectors using a quaternion"
function rotate_basis(basis,q::Quaternion)
    return [rotate(basis[1],q),
            rotate(basis[2],q),
            rotate(basis[3],q)]
end

"Returns true if within boundaries defined by maxval"
function issafe(index,maxval,minval=0)
    if index <= maxval && index > minval
        return true
    else
        return false
    end
end

"Returns true if index is within dimensions of grid"
is_safe(ind,dims) = all(map((x,y)-> x>0 && x<=y,ind,dims))

"Linear interpolation between min_val and max_val with val between 0 and 1"
linear_interp(min_val,max_val,val) = min_val + (max_val-min_val)*val

"init_val: starting value of the pixel, max_val: maximum value of the pixel, new_val: value between 0 and 1 to interpolate, cur_val: current value of the pixel"
function set_value!(init_val,max_val,new_val,cur_val)
    val = linear_interp(init_val,max_val,new_val)
    return min(cur_val+val,max_val)
end

"Linear interpolate. Returns value between 0 and 1"
function lerp(min,max,val) 
    return (val-min)/(max-min)
end

"Modified sign() that returns one if input is 0"
function dsign(val)
    if val == 0 
        return 1
    else
        return sign(val)
    end
end

"custom function to determine if val is near the target"
function isapprx(val,target,rtol)
    if (val > target*(1+rtol)) | (val < target*(1-rtol))
        return false
    else
        return true
    end
end

is_above(val,target,rtol) = (val > target) & (val < target*(1+rtol))

is_below(val,target,rtol) = (val < target) & (val > target*(1-rtol))

"find one value from a range above and one below equal to target within tolerance"
function valfound(range,target,rtol)
    return any([is_above(r,target,rtol) for r in range]) &
     any([is_below(r,target,rtol) for r in range])
end

period_avg(data,period) = sum(data[end-period:end])/(period+1)

"returns gradient of linear fit"
function lin_ext(xdata,ydata)
    if length(ydata) > 2
        coefficients  = hcat(xdata,ones(length(xdata))) \ ydata
        m = coefficients[1]
        return m
    else
        m = (ydata[2] - ydata[1]) / (xdata[2] - xdata[1])   
        return m  
    end
end

"Returns xvalue corresponding to ytarget using linear fit with gradient m"
function invert_linear(m,xpos,ypos,ytarget)
    c = ypos - m * xpos
    return (ytarget - c) / m  
end

"Take a local sample of a grid and return average value"
function sample(grid::AbstractArray,point::Vector,span::Vector)
    num = 0
    total = 0

    ndims = length(size(grid))
    ranges = Tuple(point[i]-span[i]:point[i]+span[i] for i in 1:ndims)
    
    for i in CartesianIndices(ranges)
        if checkbounds(Bool,grid,i)
            total += grid[i]
            num += 1
        end
    end
    return total / num
end

"Loop through grid and take a low resolution sample at each new grid point"
function lower_resolution(grid::AbstractArray,newdims::Tuple,span::Vector) 
    steps = cld.(size(grid),newdims)
    newgrid = zeros(eltype(grid),newdims)

    for i in CartesianIndices(newgrid)
        old_coords = (i.I .- 1) .* steps .+ 1 .+ span
        newgrid[i] = sample(grid,[old_coords...],span)
    end
    return newgrid
end

"Call lower_resolution with a constant factor to divide grid size by"
function lower_resolution(grid::AbstractArray,factor::Number) 
    newdims = cld.(size(grid),factor)
    span = [fld(factor,2) for i in 1:length(size(grid))]
    
    return lower_resolution(grid,newdims,span)
end

"Create path if it does not exist"
function make_path(path::AbstractString)
    if !isdir(path)
        mkpath(path)
    else
        error("Path $(path) already exists. Please choose a different path.")
    end
end