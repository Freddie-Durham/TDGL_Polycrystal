using TDGL_Polycrystal
using HDF5
using CUDA

function key_junction_metadata!(sim_grid, start_α, start_β, start_m⁻¹X, start_m⁻¹Y, metadata)
    sim_grid["α"] = start_α
    sim_grid["β"] = start_β
    sim_grid["m⁻¹_x"] = start_m⁻¹X
    sim_grid["m⁻¹_y"] = start_m⁻¹Y

    for (key, val) in metadata 
        HDF5.attributes(sim_grid)[key] = val
    end
end

function run_simulation(B, folder, init_time;
    pixels_per_xi, tstep, GL, levelcount, tol, conductivity, norm_resist, norm_mass,
    ramp_mode, Ecrit, Jramp, J_initial, holdtime, init_hold, grain_size, thickness,
    alphaN, betaN, init_alpha, init_beta, backend,
    rng_seed, save_states, save_frequency, resume_states, kwargs...)

    b_int = convert(Int64, round(B))
    campaign = "$(b_int)_$(b_int)%B_step1"
    filepath = "$(folder)$(campaign).h5"
    checkpoint_file = "$(folder)$(campaign)_checkpoint.h5"
    param_file = "$(folder)params$(b_int).h5"

    file_exists = isfile(filepath)
    checkpoint_exists = isfile(checkpoint_file)

    # ensure that completed simulations are never overwritten
    if file_exists && !checkpoint_exists
        @info "Simulation with magnetic field $B already complete"
        return

    # this shouldn't happen
    elseif !file_exists && checkpoint_exists
        @info "Data file missing for magnetic field $B"
        return
    end

    B /= 100 # convert from percentage to fraction of Bc2

    FindType = JcFinder
    pattern = TDGL_Polycrystal.JosephsonJunction(thickness, grain_size)

    xmin = grain_size + thickness
    ymin = grain_size + 2 * thickness
    yperiodic = false

    finder, metadata, start_α, start_β, start_m⁻¹X, start_m⁻¹Y = simulation_setup(
    pixels_per_xi, pattern, tstep, GL, conductivity, norm_resist, norm_mass, Ecrit, Jramp,
    holdtime, init_hold, xmin, ymin, yperiodic,alphaN, betaN, init_alpha, init_beta, FindType,
    levelcount, tol, backend, rng_seed, J_initial, B, ramp_mode)

    if !file_exists
        h5open(filepath, "w") do fid
            sim_grid = create_group(fid,"grid")
            key_junction_metadata!(sim_grid, start_α, start_β, start_m⁻¹X, start_m⁻¹Y, metadata)
            
            campaign_group = create_group(fid, "data")
            create_group(campaign_group, "$(B)b data")
        end
    end

    if save_states && !checkpoint_exists
        h5open(checkpoint_file, "w") do fid
            create_group(fid, "save_state")
        end
    end

    h5open(param_file, "w") do fid
        current_B = create_group(fid, "params")
        current_B["Bs"] = Vector{Float64}([B])
    end

    @info "Starting simulation with magnetic field $B" device=CUDA.device()

    finder = TDGL_Polycrystal.new_finder(finder, FindType, Ecrit, holdtime, init_hold,
             Jramp, J_initial, B, tol, levelcount, backend, rng_seed, ramp_mode)

    if resume_states && checkpoint_exists
        h5open(checkpoint_file, "r") do fid
            state_group = fid["save_state"]
            if haskey(state_group, "ψ")
                @info "Resuming from saved state for magnetic field $(B)"
                TDGL_Polycrystal.load_state!(finder, state_group)
            else
                # this also shouldn't happen
                error("Checkpoint file empty for magnetic field $(B)")
            end
        end
    end

    header = ["Current", "Electric Field", "Magnetic Field"]

    verbose_context = "B=$(B), device=$(CUDA.device())"
    sim_data, timetaken = find_jc(finder, filepath=filepath, checkpoint_file=checkpoint_file,
        save_states=save_states, save_frequency=save_frequency, B=B, header=header,
        verbose_context=verbose_context, report_period=1000)

    h5open(filepath, "r+") do fid
        campaign_group = fid["data"]
        data_group = campaign_group["$(B)b data"]
        shot_metadata = Dict("Applied field" => B, "Time taken" => timetaken)
        
        TDGL_Polycrystal.save_data(header, sim_data, data_group)
        for (key, val) in shot_metadata 
            HDF5.attributes(data_group)[key] = val
        end

        HDF5.attributes(fid["grid"])["WallTime"] = time() - init_time
        if save_states || resume_states
            rm(checkpoint_file)
        end
    end
end

function setup_simulation(;path, uID, startB, stopB, stepB, kwargs...)
    init_time = time()

    name = "$(uID)/"
    folder = path * name
    mkpath(folder)

    Bs = collect(startB:stepB:startB + stopB)
    for i in eachindex(Bs)   
        run_simulation(Bs[i], folder, init_time; kwargs...)
    end
    
    println("Simulation complete, time taken = $(time() - init_time)")
end

function main()
    kwargs = parse_CL()
    setup_simulation(;kwargs...)
end
main()