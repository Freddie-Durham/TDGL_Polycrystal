using TDGL_Polycrystal
using HDF5

function run_simulation(;uID,startB,vary_param,num_vary,
    pixels_per_xi,AA_factor,tstep,GL,levelcount,tol,conductivity,norm_resist,norm_mass,
    Ecrit,Jramp,holdtime,init_hold,N_value,rep_grain,thickness,
    xmin,ymin,yperiodic,alphaN,betaN,init_alpha,init_beta,backend,rng_seed,kwargs...)

    FindType = JcFinder
    pattern = TDGL_Polycrystal.TruncOct(N_value,xmin,rep_grain,thickness,AA_factor)

    path = "outputs/$(uID)_"*vary_param*"_converge/"
    TDGL_Polycrystal.make_path(path)
    total_time = 0

    if uppercase(vary_param) == "RNG"
        fullname = uID*"RNGseed"

        finder, metadata, start_α,start_β,start_m⁻¹,start_σ = simulation_setup(
        pixels_per_xi,pattern,tstep,GL,conductivity,norm_resist,norm_mass,Ecrit,Jramp,
        holdtime,init_hold,xmin,ymin,yperiodic,
        alphaN,betaN,init_alpha,init_beta,FindType,levelcount,tol,backend,rng_seed,
        startB)

        save_metadata(path,fullname,metadata,start_α,start_β,start_m⁻¹,start_σ)
        println("Setup Complete.")

        for _ in 1:num_vary
            println("Running simulation with rng seed = $rng_seed")

            finder = TDGL_Polycrystal.new_finder(finder,FindType,Ecrit,holdtime,init_hold,Jramp,startB,tol,levelcount,backend,rng_seed)
            sim_data, timetaken = find_jc(finder)

            filepath = "$(path)$(fullname).h5"
            h5open(filepath,"cw") do fid
                data = create_group(fid,"Seed_$rng_seed")
                data["Jc"] = sim_data[1][end]
                HDF5.attributes(data)["Bfield"] = sim_data[3][end]
                HDF5.attributes(data)["WallTime"] = timetaken
            end
            
            total_time += timetaken
            rng_seed += 1
        end
    else

        for i in 1:num_vary
            finder, metadata, start_α,start_β,start_m⁻¹,start_σ = simulation_setup(
            pixels_per_xi,pattern,tstep,GL,conductivity,norm_resist,norm_mass,Ecrit,Jramp,
            holdtime,init_hold,xmin,ymin,yperiodic,
            alphaN,betaN,init_alpha,init_beta,FindType,levelcount,tol,backend,rng_seed,
            startB) #<- last line contains arguments specific to FindType

            fullname = uID*"_"*uppercase(vary_param)*"-$i"

            save_metadata(path,fullname,metadata,start_α,start_β,start_m⁻¹,start_σ)

            println("Setup Complete.")

            sim_data, timetaken = find_jc(finder)
            header = ["Current","Electric Field","Magnetic Field"]
            save_simdata(path,fullname,sim_data,header,timetaken)
            total_time += timetaken

            if i != num_vary
                if uppercase(vary_param) == "PIXELS"
                    pixels_per_xi *= 2
                    println("Running simulation with 1/h = $pixels_per_xi")
                elseif uppercase(vary_param) == "LENGTH"
                    xmin *= 2
                    rep_grain *= 2
                    pattern = TDGL_Polycrystal.TruncOct(N_value,xmin,rep_grain,thickness,AA_factor)
                    println("Running simulation with length = $xmin")
                elseif uppercase(vary_param) == "WIDTH"
                    ymin *= 2
                    println("Running simulation with width = $ymin")
                elseif uppercase(vary_param) == "AREA"
                    xmin *= 2
                    ymin *= 2
                    rep_grain *= 2
                    pattern = TDGL_Polycrystal.TruncOct(N_value,xmin,rep_grain,thickness,AA_factor)
                    println("Running simulation with area = $(xmin*ymin)")
                elseif uppercase(vary_param) == "TOL"
                    tol *= 0.1
                    println("Running simulation with tolerance = $tol")
                elseif uppercase(vary_param) == "HOLD"
                    holdtime *= 2
                    println("Running simulation with hold time = $holdtime")
                else
                    println("Could not identify parameter: "*vary_param)
                    break
                end
            end
        end
    end

    println("Simulation complete, time taken = $(total_time)")  
end

function main()
    kwargs = parse_CL()
    run_simulation(;kwargs...)
end
main()