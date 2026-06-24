using TDGL_Polycrystal
using HDF5
using CUDA
using Base.Threads

function run_simulation(device, B, folder; kwargs...)
    CUDA.device!(device)

    @info "Starting simulation" device=CUDA.device()

    N = 10^7
    x = CUDA.fill(Float32(B), N)
    y = CUDA.rand(Float32, N)

    for _ in 1:100
        @. x = sin(x) + y
    end

    h5open("$(folder)params$(convert(Int64, round(B))).h5", "w") do fid
        current_B = create_group(fid, "params")
        current_B["Bs"] = Array(x[1:5])
    end
end

function setup_simulation(;path, uID, startB, stopB, stepB, kwargs...)
    init_time = time()
    
    name = "$(uID)/"
    folder = path * name
    mkpath(folder)

    stopB = startB + stopB
    Bs = collect(startB:stepB:stopB)
    devices = collect(CUDA.devices())

    @assert length(Bs) == length(devices) "number of simulations != number of GPUs available"

    @sync for i in eachindex(Bs)
        Threads.@spawn begin
            run_simulation(devices[i], Bs[i], folder; kwargs...)
        end
    end
    
    println("Simulation complete, time taken = $(time()-init_time)")  
end

function main()
    kwargs = parse_CL()
    setup_simulation(;kwargs...)
end
main()