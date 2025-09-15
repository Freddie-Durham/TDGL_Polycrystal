using TDGL_Polycrystal

function main()
    package_dir = abspath(joinpath(@__DIR__,".."))*"neper/2D_crystal"

    for id in 1:2
        for i in 10:500
            filename = package_dir*"D$(i)AR100xyID$id"
            if !isfile(filename*".tess")
                run(`neper -T -n from_morpho -morpho "graingrowth($(i/1000))" -dim 2 -periodicity 1 -id $id -o $filename`)
            end
        end
    end
end
main()
