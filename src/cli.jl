"parse arguments from command line, if nothing is passed in, a default value is used instead"
function parse_CL()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--path"
            help = "Specify path to store output files"
            arg_type = String
            default = "outputs/"

        "--uID"
            help = "Unique identifier, appended to data files to differentiate them"
            arg_type = String
            default = "default"

        "--startB"
            help = "Initial applied magnetic field strength"
            arg_type = Float64
            default = 0.0

        "--stopB"
            help = "Final applied magnetic field strength, expressed as a difference from initial applied field"
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

        "--anisotropy"
            help = "Anisotropy factor that increases effective mass in the z direction - used for simulating ReBCO"
            arg_type = Float64
            default = 1.0

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
            help = "Thickness of polycrystal grain boundaries in coherence lengths"
            arg_type = Float64
            default = 3.0

        "--grain_size"
            help = "Average grain size in coherence lengths"
            arg_type = Float64
            default = 32.0

        "--xmin"
            help = "Width in coherence lengths of system in the x direction"
            arg_type = Int64
            default = 64

        "--ymin"
            help = "Width in coherence lengths of system in the y direction"
            arg_type = Int64
            default = 64

        "--zmin"
            help = "Width in coherence lengths of system in the z direction"
            arg_type = Int64
            default = 64

        "--yperiodic"
            help = "Boolean value determining whether the system is periodic in the y direction (x direction, along which current is applied, is always periodic)"
            arg_type = Bool
            default = true

        "--zperiodic"
            help = "Boolean value determining whether the system is periodic in the z direction (x direction, along which current is applied, is always periodic)"
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

        "--rng_seed"
            help = "Seed for the random number generator, used to generate the initial state of the system"
            arg_type = Int64
            default = 1

        "--voronoi_seed"
            help = "Seed for the Voronoi grain distribution, used to generate the Voronoi tessellation"
            arg_type = Int64
            default = 1

        "--dims"
            help = "Number of dimensions for the simulation (2 or 3)"
            arg_type = Int64
            default = 2
        
        "--J_initial"
            help = "Initial current density to start Jc finding from (defaults to 0.0001)"
            arg_type = Float64
            default = 0.0001

        "--ramp_mode"
            help = "Mode of operation for finding Jc (default is exponential current ramp)"
            arg_type = String
            default = "EXP_INCREASE"
        
        "--save_states"
            help = "Boolean value determining whether to save the simulation state after each B field step"
            arg_type = Bool
            default = false

        "--save_frequency"
            help = "Number of time-steps between each saved state, only relevant if --save_states is true"
            arg_type = Int64
            default = 10_000

        "--resume_states"
            help = "Boolean value determining whether to resume from the last saved state or start anew"
            arg_type = Bool
            default = false

        "--continuous"
            help = "Boolean value determining whether to continue from the state reached at the previous B field or start anew"
            arg_type = Bool
            default = false
    end
    return parse_args(s, as_symbols=true)
end