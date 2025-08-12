using TDGL_Polycrystal

"run simulation setup and save metadata to file. use find_jc to run simulation"
function run_simulation(;uID,startB,stopB,stepB,num_samples,
    pixels_per_xi,AA_factor,tstep,GL,levelcount,tol,conductivity,norm_resist,norm_mass,
    Ecrit,Jramp,holdtime,init_hold,N_value,rep_grain,thickness,
    xmin,ymin,yperiodic,alphaN,betaN,init_alpha,init_beta,backend,kwargs...)
    
    FindType = BVarLinXFinder

    finder, metadata, start_α, start_β, start_m⁻¹,start_σ = simulation_setup(
    pixels_per_xi,AA_factor,N_value,rep_grain,thickness,
    tstep,GL,conductivity,norm_resist,norm_mass,
    Ecrit,Jramp,holdtime,init_hold,xmin,ymin,
    yperiodic,alphaN,betaN,init_alpha,init_beta,FindType,levelcount,
    tol,backend,
    startB,stopB,stepB,num_samples) #<- last line contains arguments specific to FindType

    path = "outputs/$(uID)LinX/"
    name = "$(uID)Bi-"*to_string(startB)*"Bf-"*to_string(stopB)
    mkpath(path*name)

    save_metadata(path,name,metadata,start_α,start_β,start_m⁻¹,start_σ)

    println("Setup Complete")

    sim_data, timetaken = find_jc(finder)

    header = ["Current","Magnetic Field","Electric Field","Current Data","E Field Data"]
    save_simdata(path,name,sim_data,header,timetaken)
    println("Simulation complete, time taken = $(timetaken)")  
end

function main()
    kwargs = parse_CL()
    run_simulation(;kwargs...)
end
main()