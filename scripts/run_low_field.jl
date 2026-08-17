using TDGL_Polycrystal
using HDF5

function run_simulation(;path, uID, startB, 
    pixels_per_xi, tstep, GL, levelcount, tol, conductivity, norm_resist, norm_mass, 
    ramp_mode, Ecrit, Jramp, J_initial, holdtime, init_hold, grain_size, thickness, 
    xmin, ymin, zmin, yperiodic, zperiodic, alphaN, betaN, init_alpha, init_beta, backend, 
    rng_seed, voronoi_seed, dims, save_states, save_frequency, resume_states, kwargs...)
    init_time = time()

    B = startB
    campaign = "e$(round(Int64, -log10(B)))"

    FindType = JcFinder
    pattern = TDGL_Polycrystal.Voronoi(grain_size, thickness, voronoi_seed)

    if dims < 3
        finder, metadata, start_α, start_β, start_m⁻¹, start_σ = simulation_setup(
        pixels_per_xi, pattern, tstep, GL, conductivity, norm_resist, norm_mass, Ecrit, Jramp,
        holdtime, init_hold, xmin, ymin, yperiodic,
        alphaN, betaN, init_alpha, init_beta, FindType, levelcount, tol, backend, rng_seed,
        J_initial, B, ramp_mode)
    else
        finder, metadata, weights = simulation_setup_3D(
        pixels_per_xi, pattern, tstep, GL, conductivity, norm_resist, norm_mass, Ecrit, Jramp, 
        holdtime, init_hold, xmin, ymin, zmin, yperiodic, zperiodic, 
        alphaN, betaN, init_alpha, init_beta, FindType, levelcount, tol, backend, rng_seed, 
        J_initial, B, ramp_mode)
    end

    name = "$(uID)/"
    mkpath(path * name)

    # save params of B field range for plotting
    h5open("$(path)$(name)params$(campaign).h5", "w") do fid
        current_B = create_group(fid, "params")
        current_B["Bs"] = [B]
    end

    filepath = "$(path)$(name)$(campaign).h5"
    checkpoint_file = "$(path)$(name)$(campaign)_checkpoint.h5"
    header = ["Current", "Electric Field", "Magnetic Field"]
    if !isfile(filepath)
        h5open(filepath, "w") do fid
            sim_grid = create_group(fid, "grid")

            # weights is either a single 3D mesh or a tuple of 2D meshes
            if dims < 3
                TDGL_Polycrystal.key_metadata(sim_grid, start_α, start_β, start_m⁻¹, start_σ, metadata)
            else    
                TDGL_Polycrystal.key_metadata(sim_grid, weights, metadata)
            end
            campaign_group = create_group(fid, "data")
            create_group(campaign_group, "$(B)b data")
        end
    end

    if save_states && !isfile(checkpoint_file)
        h5open(checkpoint_file, "w") do fid
            create_group(fid, "save_state")
        end
    end
    println("Setup Complete")

    finder = TDGL_Polycrystal.new_finder(
        finder, FindType, Ecrit, holdtime, init_hold, Jramp, J_initial, B, tol, levelcount, backend, rng_seed, ramp_mode)
    
    if finder.ramp_method isa TDGL_Polycrystal.Ramp_Method{TDGL_Polycrystal.Exp_Decrease}
        println("Ramping down current")
        J_initial = max(J_initial * 0.8, Ecrit * 1.25) # should account for conductivity
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
    sim_data, timetaken = find_jc(finder, filepath=filepath, checkpoint_file=checkpoint_file,
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