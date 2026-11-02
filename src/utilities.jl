#alias for 3D Static Vectors
const vec3 = SVector{3,Float64}
#allows broadcast dict like function
(d::Dict)(k) = d[k]

abstract type Finder end

struct JcInitHold end
struct JcJHold end
struct JcDone end

Base.string(::JcInitHold) = "Initial Hold"
Base.string(::JcJHold) = "Hold at fixed J"
Base.string(::JcDone) = "Simulation Complete"

state(finder::Finder) = MulTDGL.state(finder.solver)
system(finder::Finder) = MulTDGL.system(finder.solver)

const JcMode = Union{JcInitHold, JcJHold, JcDone}

"Linear interpolation between min_val and max_val with val between 0 and 1"
linear_interp(min_val, max_val, val) = min_val + (max_val - min_val) * val

MulTDGL_backend(d::MulTDGL.AbstractForm) = KernelAbstractions.get_backend(data(d))

const VERBOSE_LOCK = ReentrantLock()

function verbose_prefix(context)
    task_label = Base.Threads.nthreads() == 1 ? "" : "thread $(Base.Threads.threadid())"

    if isempty(context)
        return isempty(task_label) ? "" : "[$task_label] "
    elseif isempty(task_label)
        return "[$context] "
    else
        return "[$context, $task_label] "
    end
end

"print a complete verbose update without interleaving concurrent threaded output"
function verbose_update(lines; context="")
    lock(VERBOSE_LOCK)
    try
        prefix = verbose_prefix(context)
        for line in lines
            println(prefix, line)
        end
        flush(stdout)
    finally
        unlock(VERBOSE_LOCK)
    end
end

"save initial grid settings and metadata to HDF5 file"
function key_metadata(sim_grid, start_α, start_β, start_m⁻¹, start_σ, metadata)
    # save initial grid settings, can save on every timestep?
    sim_grid["α"] = start_α
    sim_grid["β"] = start_β
    sim_grid["m⁻¹_x"] = start_m⁻¹
    sim_grid["σ_y"] = start_σ

    # append metadata for whole campaign to grid
    for (key, val) in metadata 
        HDF5.attributes(sim_grid)[key] = val
    end
end

"save initial grid settings (only weights) and metadata to HDF5 file"
function key_metadata(sim_grid, weights, metadata)
    # save initial grid settings, can save on every timestep?
    sim_grid["Weights"] = weights

    # append metadata for whole campaign to grid
    for (key, val) in metadata 
        HDF5.attributes(sim_grid)[key] = val
    end
end

"create path to HDF5 file and save metadata and initial grid settings"
function save_metadata(path, name, metadata, start_α, start_β, start_m⁻¹, start_σ)
    filepath = "$(path)$(name).h5"
    h5open(filepath, "cw") do fid
        sim_grid = create_group(fid, "grid")
        key_metadata(sim_grid, start_α, start_β, start_m⁻¹, start_σ, metadata)
    end
end

function save_state(finder::Finder, state_group, current, B)
    s = state(finder)
    for (lbl, vals) in zip(["ψ", "a", "φ"], [s.ψ, s.a, s.φ])
        if haskey(state_group, lbl)
            delete_object(state_group, lbl)
        end
        state_group[lbl] = adapt(CPU(), data(vals))
    end

    if haskey(HDF5.attributes(state_group), "Current")
        delete!(attrs(state_group), "Current")
    end
    HDF5.attributes(state_group)["Current"] = current

    if haskey(HDF5.attributes(state_group), "Magnetic Field")
        delete!(attrs(state_group), "Magnetic Field")
    end
    HDF5.attributes(state_group)["Magnetic Field"] = B

    # progress of the finder, so a resumed run continues exactly where it stopped
    # rather than restarting the initial hold / current hold
    progress = Dict{String,Any}()
    hasproperty(finder, :mode) && (progress["Mode"] = string(finder.mode))
    hasproperty(finder, :curholdsteps) && (progress["Hold steps"] = finder.curholdsteps)
    hasproperty(finder, :esum) && (progress["E sum"] = finder.esum)
    # find_jc records E_field *before* each step, so without this the first row
    # after a resume would be the freshly-initialised value 0
    hasproperty(finder, :E_field) && (progress["Electric Field"] = finder.E_field)
    for (key, val) in progress
        haskey(HDF5.attributes(state_group), key) && delete!(attrs(state_group), key)
        HDF5.attributes(state_group)[key] = val
    end
end

"finder modes that can be restored from a checkpoint, keyed by `string(mode)`"
const RESTORABLE_MODES = Dict(string(m) => m for m in (JcInitHold(), JcJHold(), JcDone()))

"loads raw data from HDF5 file into the state of the finder"
function load_state!(finder::Finder, state_group)
    s = state(finder)
    backend = MulTDGL_backend(s.ψ)
    for lbl in ["ψ", "a", "φ"]
        if haskey(state_group, lbl)
            target = getproperty(s, Symbol(lbl))
            data_array = read(state_group[lbl])

            data(target) .= adapt(backend, data_array)
        else
            error("State group does not contain $lbl")
        end
    end

    if haskey(HDF5.attributes(state_group), "Current")
        finder.j = HDF5.read_attribute(state_group, "Current")
    end

    if haskey(HDF5.attributes(state_group), "Magnetic Field")
        finder.B_field = HDF5.read_attribute(state_group, "Magnetic Field")
        # the field is applied through the system, so re-apply it if it was being ramped
        if hasproperty(finder, :ramp_method) && finder.ramp_method.mode isa Ramp_B
            apply_B_field!(finder, finder.B_field)
        end
    end

    # restore the finder's progress (absent in checkpoints written by older code,
    # in which case the run restarts with the initial hold as before)
    state_attrs = HDF5.attributes(state_group)
    if haskey(state_attrs, "Mode") && hasproperty(finder, :mode)
        mode = HDF5.read_attribute(state_group, "Mode")
        if haskey(RESTORABLE_MODES, mode)
            finder.mode = RESTORABLE_MODES[mode]
        else
            @warn "Unknown mode \"$mode\" in checkpoint; restarting from the initial hold"
        end
    end
    if haskey(state_attrs, "Hold steps") && hasproperty(finder, :curholdsteps)
        finder.curholdsteps = HDF5.read_attribute(state_group, "Hold steps")
    end
    if haskey(state_attrs, "E sum") && hasproperty(finder, :esum)
        finder.esum = HDF5.read_attribute(state_group, "E sum")
    end
    if haskey(state_attrs, "Electric Field") && hasproperty(finder, :E_field)
        finder.E_field = HDF5.read_attribute(state_group, "Electric Field")
    end

    # E computed from the change in A since the previous step (SUBSECTION = true in
    # FindJc.jl) needs the previous A, which is the A just loaded. Only copied when
    # needed, as a full state copy is large on the GPU.
    if SUBSECTION && hasproperty(finder, :prev_state)
        finder.prev_state = copy(s)
    end
end

"iterate through header and save data to HDF5 file"
function save_data(header, sim_data, data_group)
    for (i, h) in enumerate(header)
        if haskey(data_group, h)
            old_data = read(data_group[h])
            append!(old_data, sim_data[i])
            delete_object(data_group, h)
            data_group[h] = old_data

        else
            data_group[h] = sim_data[i]
        end
    end
end

"save results of simulation and time taken to same HDF5 file"
function save_simdata(path, name, sim_data, header, walltime)
    filepath = "$(path)$(name).h5"
    h5open(filepath,"cw") do fid
        data_group = create_group(fid,"data")
        save_data(header, sim_data, data_group)
        sim_grid = fid["grid"]
        HDF5.attributes(sim_grid)["Wall Time"] = walltime
    end
end

"""
    atomic_h5_update(f, path; copy_existing=true)

Crash-safe update of the HDF5 file at `path`. `f(fid)` is applied to a scratch file
`path * ".tmp"` (a copy of `path` if `copy_existing`, otherwise a new empty file),
which is flushed to disk and then atomically renamed over `path`. If the process is
killed at any point, `path` holds either the complete old or the complete new
version, never a half-written one. Returns the value of `f(fid)`.
"""
function atomic_h5_update(f, path; copy_existing=true)
    tmp = path * ".tmp"
    if copy_existing && isfile(path)
        cp(path, tmp; force=true)
        result = h5open(f, tmp, "r+")
    else
        result = h5open(f, tmp, "w")
    end

    # make sure the data has reached the disk before the rename makes it visible
    open(tmp, "r") do io
        rawfd = fd(io)  # an Int on Julia 1.11, a RawFD on newer versions
        cfd = rawfd isa Integer ? Cint(rawfd) : reinterpret(Cint, rawfd)
        ccall(:fsync, Cint, (Cint,), cfd) == 0 || @warn "fsync failed for $tmp"
    end
    # rename(2) atomically replaces `path` (unlike `mv(...; force=true)`, which
    # deletes the destination first)
    Base.Filesystem.rename(tmp, path)
    return result
end

"""
    write_checkpoint(finder, checkpoint_file, data_rows)

Atomically replace `checkpoint_file` with the current state of `finder`, recording
`data_rows`, the number of rows the data file holds at this point of the simulation.
"""
function write_checkpoint(finder::Finder, checkpoint_file, data_rows)
    atomic_h5_update(checkpoint_file; copy_existing=false) do fid
        state_group = create_group(fid, "save_state")
        save_state(finder, state_group, finder.j, finder.B_field)
        HDF5.attributes(state_group)["Data rows"] = data_rows
    end
end

"number of rows saved under `header` in `data_group` (0 if nothing saved yet)"
function saved_rows(header, data_group)
    (isempty(header) || !haskey(data_group, header[1])) && return 0
    return length(data_group[header[1]])
end

"""
    sync_data_with_checkpoint(filepath, checkpoint_file, data_path, header)

Called when resuming. The data file is updated before the checkpoint, so a job killed
between the two leaves rows in the data file that the checkpoint state does not know
about; these would be duplicated when the simulation resumes. Truncate the data under
`data_path` back to the row count stored in the checkpoint.
"""
function sync_data_with_checkpoint(filepath, checkpoint_file, data_path, header)
    (isempty(header) || !isfile(checkpoint_file) || !isfile(filepath)) && return

    expected = h5open(checkpoint_file, "r") do fid
        state_group = fid["save_state"]
        haskey(HDF5.attributes(state_group), "Data rows") || return nothing
        HDF5.read_attribute(state_group, "Data rows")
    end
    isnothing(expected) && return  # checkpoint from older code, or no state saved yet

    actual = h5open(fid -> saved_rows(header, fid[data_path]), filepath, "r")
    if actual > expected
        atomic_h5_update(filepath) do fid
            data_group = fid[data_path]
            for h in header
                haskey(data_group, h) || continue
                vals = read(data_group[h])[1:expected]
                delete_object(data_group, h)
                data_group[h] = vals
            end
        end
    elseif actual < expected
        @warn "Data file has fewer rows ($actual) than the checkpoint expects ($expected)"
    end
end
