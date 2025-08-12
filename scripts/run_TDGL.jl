using TDGL_Polycrystal
using HDF5

function run_simulation(;uID,startB,stopB,stepB,
    pixels_per_xi,AA_factor,tstep,GL,levelcount,tol,conductivity,norm_resist,norm_mass,
    Ecrit,Jramp,holdtime,init_hold,N_value,rep_grain,thickness,
    xmin,ymin,yperiodic,alphaN,betaN,init_alpha,init_beta,backend,kwargs...)
    init_time = time()

    stopB = startB + stopB

    #avoid difficulty with filenames that have decimal points
    if stepB >= 1 
        campaign = "$(convert(Int64,startB))_$(convert(Int64,stopB))%B_step$(convert(Int64,stepB))"
    else
        inv_step = convert(Int64,round(1/stepB))
        campaign = "$(convert(Int64,round(startB)))_$(convert(Int64,round(stopB)))%B_step1_$(inv_step)"
    end

    B_range = (startB/100):(stepB/100):(stopB/100)

    FindType = JC2DFinder
    finder, metadata, start_α,start_β,start_m⁻¹,start_σ = simulation_setup(
    pixels_per_xi,AA_factor,N_value,rep_grain,thickness,
    tstep,GL,conductivity,norm_resist,norm_mass,
    Ecrit,Jramp,holdtime,init_hold,xmin,ymin,
    yperiodic,alphaN,betaN,init_alpha,init_beta,FindType,levelcount,tol,backend,
    B_range[1])

    path = "outputs/"
    name = "$(uID)/"
    mkpath(path*name)

    #Save params of B field range for plotting
    h5open("$(path)$(name)params$(convert(Int64,round(startB))).h5","w") do fid
        current_B = create_group(fid,"params")
        current_B["Bs"] = Vector{Float64}(B_range)
    end

    filepath = "$(path)$(name)$(campaign).h5"
    header = ["Current","Electric Field"]
    h5open(filepath,"w") do fid
        sim_grid = create_group(fid,"grid")
        TDGL_Polycrystal.key_metadata(sim_grid,start_α,start_β,start_m⁻¹,start_σ,metadata)
        println("Setup Complete")

        #create folder for current campaign
        campaign_group = create_group(fid,"data")
    
        #iterate through B fields, recording data and shot-specific metadata
        for B in B_range
            finder, α, β, m⁻¹,σ = simulation_setup(
            pixels_per_xi,AA_factor,N_value,rep_grain,thickness,
            tstep,GL,conductivity,norm_resist,norm_mass,
            Ecrit,Jramp,holdtime,init_hold,xmin,ymin,
            yperiodic,alphaN,betaN,init_alpha,init_beta,FindType,levelcount,tol,backend,
            B)

            println("Running simulation with B = $(B)")
            sim_data, timetaken = find_jc(finder)

            shot_metadata = Dict("Applied field" => B,"Time taken" => timetaken)
            data_group = create_group(campaign_group,"$(B)b data")
            TDGL_Polycrystal.save_data(header,sim_data,data_group)
            for (key,val) in shot_metadata 
                HDF5.attributes(data_group)[key] = val
            end
        end 
        HDF5.attributes(sim_grid)["WallTime"] = time()-init_time
    end 
    println("Simulation complete, time taken = $(time()-init_time)")  
end

function main()
    kwargs = parse_CL()
    run_simulation(;kwargs...)
end
main()

