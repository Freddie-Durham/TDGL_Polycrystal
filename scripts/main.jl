using TDGL_Polycrystal

function generate_2D_tess(dir)
    dir *= "2D_crystal/"

    for id in 1:2
        for i in 5:10
            filename = dir*"D$(i)AR100xyID$id"
            if !isfile(filename*".tess")
                run(`neper -T -n from_morpho -morpho "graingrowth($(i/1000))" -dim 2 -periodicity 1 -id $id -o $filename`)
            end
        end
    end
end

function generate_3D_tess(dir)
    dir *= "3D_crystal/"
    xdims = 150

    for seed in 1:2
        for D in [22.4,28.8,35.2,41.6,48.0,54.4]
            grain_size = D/xdims
            TDGL_Polycrystal.create_3D_tessellation(dir,grain_size,seed)
        end
    end
end

function main()
    dir = abspath(joinpath(@__DIR__,".."))*"neper/"
    generate_3D_tess(dir)
end
main()
