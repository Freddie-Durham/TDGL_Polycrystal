using TDGL_Polycrystal
using HDF5

function run_simulation(;path, uID, startB, stopB, stepB, 
    pixels_per_xi, tstep, GL, levelcount, tol, conductivity, norm_resist, norm_mass, 
    ramp_mode, Ecrit, Jramp, J_initial, holdtime, init_hold, grain_size, thickness, 
    xmin, ymin, zmin, yperiodic, zperiodic, alphaN, betaN, init_alpha, init_beta, backend, 
    rng_seed, voronoi_seed, dims, save_states, save_frequency, resume_states, continuous, kwargs...)
    init_time = time()

    stopB = startB + stopB

    #avoid difficulty with filenames that have decimal points
    if stepB >= 1 
        campaign = "$(convert(Int64,startB))_$(convert(Int64, stopB))%B_step$(convert(Int64, stepB))"
    else
        inv_step = convert(Int64, round(1 / stepB))
        campaign = "$(convert(Int64, round(startB)))_$(convert(Int64, round(stopB)))%B_step1_$(inv_step)"
    end

    B_range = (startB / 100):(stepB / 100):(stopB / 100)

    FindType = JcFinder
    pattern = TDGL_Polycrystal.Voronoi(grain_size, thickness, voronoi_seed)

    if dims < 3
        finder, metadata, start_α, start_β, start_m⁻¹, start_σ = simulation_setup(
        pixels_per_xi, pattern, tstep, GL, conductivity, norm_resist, norm_mass, Ecrit, Jramp,
        holdtime, init_hold, xmin, ymin, yperiodic,
        alphaN, betaN, init_alpha, init_beta, FindType, levelcount, tol, backend, rng_seed,
        J_initial, B_range[1], ramp_mode)
    else
        finder, metadata, weights = simulation_setup_3D(
        pixels_per_xi, pattern, tstep, GL, conductivity, norm_resist, norm_mass, Ecrit, Jramp, 
        holdtime, init_hold, xmin, ymin, zmin, yperiodic, zperiodic, 
        alphaN, betaN, init_alpha, init_beta, FindType, levelcount, tol, backend, rng_seed, 
        J_initial, B_range[1], ramp_mode)
    end

    name = "$(uID)/"
    mkpath(path*name)

    #Save params of B field range for plotting
    h5open("$(path)$(name)params$(convert(Int64,round(startB))).h5", "w") do fid
        current_B = create_group(fid, "params")
        current_B["Bs"] = Vector{Float64}(B_range)
    end

    filepath = "$(path)$(name)$(campaign).h5"
    checkpoint_file = "$(path)$(name)$(campaign)_checkpoint.h5"
    header = ["Current", "Electric Field", "Magnetic Field"]
    if !isfile(filepath)
        h5open(filepath, "w") do fid
            sim_grid = create_group(fid,"grid")

            #weights is either a single 3D mesh or a tuple of 2D meshes
            if dims < 3
                TDGL_Polycrystal.key_metadata(sim_grid, start_α, start_β, start_m⁻¹, start_σ, metadata)
            else    
                TDGL_Polycrystal.key_metadata(sim_grid, weights, metadata)
            end
            campaign_group = create_group(fid, "data")
            for B in B_range
                create_group(campaign_group, "$(B)b data")
            end
        end
    end

    if save_states && !isfile(checkpoint_file)
        h5open(checkpoint_file, "w") do fid
            create_group(fid, "save_state")
        end
    end
    println("Setup Complete")

    #iterate through B fields, recording data and shot-specific metadata
    for B in B_range
        if !continuous || B == B_range[1]
            finder = TDGL_Polycrystal.new_finder(
            finder, FindType, Ecrit, holdtime, init_hold, Jramp, J_initial, B, tol, levelcount, backend, rng_seed, ramp_mode)
            
            if finder.ramp_method isa TDGL_Polycrystal.Ramp_Method{TDGL_Polycrystal.Exp_Decrease}
                println("Ramping down current")
                J_initial = max(J_initial * 0.8, Ecrit * 1.25) #should account for conductivity
            end
        else
            apply_B_field!(finder, B)
            finder.j *= 1.2
            finder.mode = JcInitHold()
        end

        if resume_states
            if isfile(checkpoint_file)
                h5open(checkpoint_file, "r") do fid
                    state_group = fid["save_state"]
                    if haskey(state_group, "ψ")
                        println("Resuming from saved state for B = $(B)")
                        TDGL_Polycrystal.load_state!(finder, state_group)
                    else
                        println("Saved state not found for B = $(B), starting new simulation")
                    end
                end
             
            else
                println("Save file not found for B = $(B), starting new simulation")   
            end
        end

        println("Running simulation with B = $(B)")
        sim_data, timetaken = find_jc(finder, filepath=filepath, checkpoint_file = checkpoint_file,
            save_states=save_states, save_frequency=save_frequency, B=B, header=header)

        h5open(filepath, "r+") do fid
            campaign_group = fid["data"]
            data_group = campaign_group["$(B)b data"]
            shot_metadata = Dict("Applied field" => B, "Time taken" => timetaken)
            
            TDGL_Polycrystal.save_data(header, sim_data, data_group)
            for (key, val) in shot_metadata 
                HDF5.attributes(data_group)[key] = val
            end
        end
    end

    h5open(filepath,"r+") do fid
        HDF5.attributes(fid["grid"])["WallTime"] = time() - init_time
        if save_states
            rm(checkpoint_file)
        end
    end
    
    println("Simulation complete, time taken = $(time()-init_time)")  
end

function main()
    kwargs = parse_CL()
    run_simulation(;kwargs...)
end
main()

