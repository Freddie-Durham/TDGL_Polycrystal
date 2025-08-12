using TDGL_Polycrystal

function run_simulation(;uID,startB,vary_param,num_vary,
    pixels_per_xi,AA_factor,tstep,GL,levelcount,tol,conductivity,norm_resist,norm_mass,
    Ecrit,Jramp,holdtime,init_hold,N_value,rep_grain,thickness,
    xmin,ymin,yperiodic,alphaN,betaN,init_alpha,init_beta,backend,kwargs...)
    
    FindType = JC2DFinder

    path = "outputs/$(uID)_"*vary_param*"_converge/"
    name = "$(uID)B-"*to_string(startB)
    TDGL_Polycrystal.make_path(path)
    total_time = 0

    for i in 1:num_vary
        finder, metadata, start_α, start_β, start_m⁻¹,start_σ = simulation_setup(
        pixels_per_xi,AA_factor,N_value,rep_grain,thickness,
        tstep,GL,conductivity,norm_resist,norm_mass,
        Ecrit,Jramp,holdtime,init_hold,xmin,ymin,
        yperiodic,alphaN,betaN,init_alpha,init_beta,FindType,levelcount,
        tol,backend,
        startB) #<- last line contains arguments specific to FindType

        fullname = name*"_"*uppercase(vary_param)*"-$i"

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
                println("Running simulation with length = $xmin")
            elseif uppercase(vary_param) == "WIDTH"
                ymin *= 2
                println("Running simulation with width = $ymin")
            elseif uppercase(vary_param) == "AREA"
                xmin *= 2
                ymin *= 2
                rep_grain *= 2
                println("Running simulation with area = $(xmin*ymin)")
            elseif uppercase(vary_param) == "TOL"
                tol *= 0.1
                println("Running simulation with tolerance = $tol")
            else
                println("Could not identify parameter: "*vary_param)
                break
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