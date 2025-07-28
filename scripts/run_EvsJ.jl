using TDGL_Polycrystal

function run_simulation(;uID,startB,max_steps,
    pixels_per_xi,AA_factor,tstep,GL,levelcount,tol,conductivity,norm_resist,norm_mass,
    Ecrit,Jramp,holdtime,init_hold,N_value,N_crystal,thickness,
    xmin,ymin,yperiodic,alphaN,betaN,backend,kwargs...)
    
    FindType = EvsJ

    finder, metadata, start_α, start_β, start_m⁻¹,start_σ = simulation_setup(
    pixels_per_xi,AA_factor,N_value,N_crystal,thickness,
    tstep,GL,conductivity,norm_resist,norm_mass,
    Ecrit,Jramp,holdtime,init_hold,xmin,ymin,
    yperiodic,alphaN,betaN,FindType,levelcount,
    tol,backend,
    startB,max_steps) #<- last line contains arguments specific to FindType

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

function main()
    kwargs = parse_CL()
    run_simulation(;kwargs...)
end
main()