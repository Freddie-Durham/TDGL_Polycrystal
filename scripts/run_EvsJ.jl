using TDGL_Polycrystal
using Profile

function run_simulation(;uID,startB,max_steps,
    pixels_per_xi,AA_factor,tstep,GL,levelcount,tol,conductivity,norm_resist,norm_mass,
    Ecrit,Jramp,holdtime,init_hold,N_value,rep_grain,thickness,
    xmin,ymin,yperiodic,alphaN,betaN,init_alpha,init_beta,backend,rng_seed,profil,kwargs...)

    FindType = EvsJ

    finder, metadata, start_α, start_β, start_m⁻¹,start_σ = simulation_setup(
    pixels_per_xi,AA_factor,N_value,rep_grain,thickness,
    tstep,GL,conductivity,norm_resist,norm_mass,
    Ecrit,Jramp,holdtime,init_hold,xmin,ymin,
    yperiodic,alphaN,betaN,init_alpha,init_beta,FindType,levelcount,
    tol,backend,rng_seed,
    startB,max_steps) #<- last line contains arguments specific to FindType

    if profil
        prof_sim(finder)
    else
        path = "outputs/$(uID)_EvsJ/"
        name = "$(uID)B-"*to_string(startB)
        TDGL_Polycrystal.make_path(path)

        save_metadata(path,name,metadata,start_α,start_β,start_m⁻¹,start_σ)

        println("Setup Complete.")

        sim_data, timetaken = find_jc(finder)
        header = ["Current","Super Current","Magnetic Field","Electric Field"]
        save_simdata(path,name,sim_data,header,timetaken)
        println("Simulation complete, time taken = $(timetaken)")  
    end
end

function prof_sim(finder)
    ###profile the simulation
    run_small(finder) 
    @profile run_small(finder)

    Profile.print()
end

function run_small(finder)
    for _ in 1:100
        step!(finder)
    end
end

function main()
    kwargs = parse_CL()
    run_simulation(;kwargs...)
end
main()