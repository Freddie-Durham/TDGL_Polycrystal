"save initial grid settings and metadata to HDF5 file"
function key_metadata(sim_grid,start_α,start_β,start_m⁻¹,start_σ,metadata)
    #save initial grid settings, can save on every timestep?
    sim_grid["α"] = start_α
    sim_grid["β"] = start_β
    sim_grid["m⁻¹"] = start_m⁻¹
    sim_grid["σ"] = start_σ

    #append metadata for whole campaign to grid
    for (key,val) in metadata 
        HDF5.attributes(sim_grid)[key] = val
    end
end

"create path to HDF5 file and save metadata and initial grid settings"
function save_metadata(path,name,metadata,start_α,start_β,start_m⁻¹,start_σ)
    filepath = "$(path)$(name).h5"
    h5open(filepath,"cw") do fid
        sim_grid = create_group(fid,"grid")
        key_metadata(sim_grid,start_α,start_β,start_m⁻¹,start_σ,metadata)
    end
end

"iterate through header and save data to HDF5 file"
function save_data(header,sim_data,data_group)
    for (i,h) in enumerate(header)
        data_group["$h"] = sim_data[i]
    end
end

"save results of simulation and time taken to same HDF5 file"
function save_simdata(path,name,sim_data,header,walltime)
    filepath = "$(path)$(name).h5"
    h5open(filepath,"cw") do fid
        data_group = create_group(fid,"data")
        save_data(header,sim_data,data_group)
        sim_grid = fid["grid"]
        HDF5.attributes(sim_grid)["Wall Time"] = walltime
    end
end

"convert number fraction to string percentage"
function to_string(num)
    num = convert(Int,floor(num*100))
    return string(num)
end

"parse arguments from command line, if nothing is passed in, a default value is used instead"
function parse_CL()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--uID"
            help = "Unique identifier, appended to data files to differentiate them"
            arg_type = String
            default = "default"

        "--startB"
            help = "Initial applied magnetic field strength"
            arg_type = Float64
            default = 0.0

        "--stopB"
            help = "Final applied magnetic field strength"
            arg_type = Float64
            default = 0.1

        "--stepB"
            help = "Magnitude by which magnetic field strength is incremented"
            arg_type = Float64
            default = 0.05

        "--num_samples"
            help = "Maximum number of E(J) measurements taken for a given applied magnetic field"
            arg_type = Int64
            default = 5

        "--max_steps"
            help = "Maximum number of steps taken in fixed J/B mode, after which the simulation is stopped"
            arg_type = Int64
            default = 1000

        "--pixels_per_xi"
            help = "Number of grid steps in a single coherence length. Grid step h=1/pixels_per_xi"
            arg_type = Int64
            default = 4

        "--AA_factor"
            help = "Factor by which resolution of grid is increased to anti-alias the rasterised pattern"
            arg_type = Int64 
            default = 5

        "--tstep"
            help = "Time increment per simulation step, in the normalised TDGL time unit τ"
            arg_type = Float64
            default = 0.5

        "--GL"
            help = "Ginzburg-Landau parameter for the superconductor"
            arg_type = Float64
            default = 10.0

        "--levelcount"
            help = "Number of discretised grids of different coarseness used in the multigrid method"
            arg_type = Int64
            default = 3

        "--tol"
            help = "Tolerance of iterative matrix method, every value in residual vector must be below tol for matrix equation to be considered solved"
            arg_type = Float64
            default = 1e-3

        "--conductivity"
            help = "Normal state conductivity of the superconductor"
            arg_type = Float64
            default = 1.0

        "--norm_resist"
            help = "Resistivity of the normal material in the superconducting matrix (proportional to effective carrier mass)"
            arg_type = Float64
            default = 2.0

        "--norm_mass"
            help = "Effective carrier mass of the normal material in the superconducting matrix (proportional to resistivity)"
            arg_type = Float64
            default = 2.0

        "--Ecrit"
            help = "Critical electric field criterion, used to determine the critical current density"
            arg_type = Float64
            default = 1e-4
        
        "--Jramp"
            help = "Magnitude by which critical current density is incremented"
            arg_type = Float64
            default = 1e-4

        "--holdtime"
            help = "Number of timesteps (--tstep) the system is held for at the same magnetic field and current density to reach equilibrium"
            arg_type = Int64
            default = 500

        "--init_hold"
            help = "Number of timesteps (--tstep) the system is held for initially to allow the condensation of the order parameter"
            arg_type = Int64
            default = 100

        "--N_value"
            help = "Determine angle θ of crystal grains with respect to the x axis using the formula θ = tanh(1/N_value)"
            arg_type = Int64
            default = 2

        "--rep_grain"
            help = "Number of times grain structure is repeated in simulation domain. Must be a factor of --xmin"
            arg_type = Int64
            default = 0

        "--thickness"
            help = "Thickness of polycrystal grain boundaries"
            arg_type = Float64
            default = 3.0

        "--xmin"
            help = "Width in coherence lengths of system in the x direction"
            arg_type = Int64
            default = 64

        "--ymin"
            help = "Width in coherence lengths of system in the y direction"
            arg_type = Int64
            default = 64

        "--yperiodic"
            help = "Boolean value determining whether the system is periodic in the y direction (x direction, along which current is applied, is always periodic)"
            arg_type = Bool
            default = true

        "--alphaN"
            help = "TDGL condensation parameter for the normal state material"
            arg_type = Float64
            default = -1.0

        "--betaN"
            help = "TDGL nonlinearity parameter for the normal state material"
            arg_type = Float64
            default = 1.0     
            
        "--init_alpha"
            help = "TDGL condensation parameter for the superconducting material"
            arg_type = Float64
            default = 1.0

        "--init_beta"
            help = "TDGL nonlinearity parameter for the superconducting material"
            arg_type = Float64
            default = 1.0 
        
        "--backend"
            help = "Type of hardware the simulation will be run on (eg. CPU, CUDA etc.)"
            arg_type = String
            default = "CPU"

        "--profil"
            help = "Debug setting, set to true to run a profiler on an EvsJ trace"
            arg_type = Bool
            default = false
        
        "--vary_param"
            help = "Debug setting, choose parameter to vary to test convergence when finding Jc"
            arg_type = String
            default = "pixels"

        "--num_vary"
            help = "Debug setting, determine number of simulations to run when varying --vary_param"
            arg_type = Int64
            default = 1
    end
    return parse_args(s,as_symbols=true)
end