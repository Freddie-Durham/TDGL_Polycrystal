using GLMakie
const rad2deg = 180/π

function test_tesselate()
    pixels = 4
    xmin = 384 
    xdims = xmin*pixels
    ydims = xmin*pixels  

 
    mesh = ones(xdims*ydims)
    mesh = reshape(mesh,xdims,ydims)

    obs = Observable(mesh)
    θ,a = periodic_crystal(2,xmin)

    tesselate!(octagon!,mesh,a*pixels,1*pixels,θ,xdims,ydims)

    repx = 1
    repy = 1

    fig = Figure()
    a = lift(m -> repeat(m, repx, repy), obs)
    ax = Axis(fig[1, 1], aspect=DataAspect(),
              limits=(1, repx*xdims, 1, repy*ydims))
    hidedecorations!(ax)
    heatmap!(ax,a,colormap=:viridis)
    save("quadGBthick.png",fig)
    display(fig)
end

function get_mesh(grid)
    colpos = Vector{Point3{Int64}}([])
    colvals = Vector{Float64}([])
    
    for c in CartesianIndices(grid)
        if grid[c] > 0
            push!(colpos,Point3{Int64}(c[1],c[2],c[3]))
            push!(colvals,grid[c])
        end
    end
    return colpos,colvals
end

function show_mesh(mesh,xdims,ydims,zdims)
    colpos,colvals = get_mesh(mesh)

    #colvals = 10*(log.(colvals))/findmax(log.(colvals))[1]

    aspect=(xdims,ydims,zdims)
    perspectiveness=0.2
    fig = Figure(; size=(800,800))
    ax1 = Axis3(fig[1, 1]; aspect, perspectiveness)
    limits!(ax1,0,xdims,0,ydims,0,zdims)
    ms = meshscatter!(ax1, colpos; markersize=0.1,
    marker=Rect3f(Vec3f(0), Vec3f(1)),color = colvals,
    colormap = :tab10, colorrange = (1, 10),
    transparency=false)
    fig
end

function slice3Dmesh(mesh)
    fig = Figure(; size=(1200,1200))
    ax1 = Axis(fig[1,1])
    ax2 = Axis(fig[1,2])
    ax3 = Axis(fig[2,1])
    ax4 = Axis(fig[2,2])
    heatmap!(ax1,mesh[:,5,:])
    heatmap!(ax2,mesh[5,:,:])
    heatmap!(ax3,mesh[:,:,5])
    heatmap!(ax4,mesh[:,11,:])
    display(fig)
end

function display_from_file(name,var)
    mesh = load(name)[var]
    xdims = size(mesh)[1]
    ydims = size(mesh)[2]
    zdims = size(mesh)[3]
    show_mesh(mesh,xdims,ydims,zdims)
end

function test_periodic(mesh)
    fig = Figure()

    xf1 = mesh[1,:,:]
    xf2 = mesh[end,:,:]
    ax1 = Axis(fig[1,1],xlabel="X faces",aspect=AxisAspect(1))
    heatmap!(ax1,xf1.+2*xf2)

    yf1 = mesh[:,1,:]
    yf2 = mesh[:,end,:]
    ax2 = Axis(fig[1,2],xlabel="Y faces",aspect=AxisAspect(1))
    heatmap!(ax2,yf1.+2*yf2)

    zf1 = mesh[:,:,1]
    zf2 = mesh[:,:,end]
    ax3 = Axis(fig[1,3],xlabel="Z faces",aspect=AxisAspect(1))
    heatmap!(ax3,zf1.+2*zf2)

    fig
end

function test_cube_mesh()
    min = 3 
    pixels = 64

    xdims = min*pixels
    zdims = min*pixels
    ydims = min*pixels

    mesh = zeros(xdims*ydims*zdims)
    mesh = reshape(mesh,xdims,ydims,zdims)

    θ = 1.910810018469853
    q = attitude_quaternion(θ,[1,1,0])

    tesselateCube!(mesh,[xdims,ydims,zdims],q,64,4,10)

    test_periodic(mesh)
    #show_mesh(mesh,xdims,ydims,zdims)
end

function test_octahedra_mesh()
    min = 6 
    pixels = 4
    width = 12

    xdims = min*pixels*width
    zdims = min*pixels*width
    ydims = min*pixels*width

    mesh = zeros(xdims*ydims*zdims)
    mesh = reshape(mesh,xdims,ydims,zdims)

    q = attitude_quaternion(1.910810018469853,[1,1,0])
    #q = attitude_quaternion(0,[1,0,0])

    #trunc_oct!(mesh,[30,30,30]*pixels,orthonorm_basis(),36*pixels,1,10)
    tesselateOct!(mesh,[xdims,ydims,zdims],q,width*pixels,1,10)

    #test_periodic(mesh)
    show_mesh(mesh,xdims,ydims,zdims)
end

function test_quaternion()
    v = [1.0,0.0,0.0]
    u = [0.0,0.0,1.0]
    q = attitude_quaternion(π/2,v)
    print(rotate(u,q))
end

function main()
    #test_tesselate()
    #test_quaternion()
    test_octahedra_mesh()
    #test_cube_mesh()
end
main()