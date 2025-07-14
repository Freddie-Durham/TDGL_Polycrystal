using HDF5
using TDGL2D

function version()
    return "0.3.1"
end

function run_simulation(uniqueID,startB,stopB,stepB,vortex_radius,tstep,GL,levelcount,tol,init_σ,norm_resist,norm_inv_mass,Ecrit,Jramp,holdtime,N,num_crystal,thickness,xmin,ymin,yperiodic,alphaN,betaN,backend)
    init_time = time()

    #avoid difficulty with filenames that have decimal points
    if stepB >= 1 
        campaign = "$(convert(Int64,startB))_$(convert(Int64,stopB))%B_step$(convert(Int64,stepB))"
    else
        inv_step = convert(Int64,round(1/stepB))
        campaign = "$(convert(Int64,round(startB)))_$(convert(Int64,round(stopB)))%B_step1_$(inv_step)"
    end

    B_range = (startB/100):(stepB/100):(stopB/100)

    FindType = JC2DFinder
    init_hold = 100
    finder, metadata, start_α,start_β,start_m⁻¹,start_σ = simulation_setup(vortex_radius,N,num_crystal,thickness,tstep,GL,
    init_σ,norm_resist,norm_inv_mass,Ecrit,Jramp,holdtime,init_hold,xmin,ymin,yperiodic,alphaN,betaN,FindType,levelcount,tol,backend,string(pkgversion(TDGL2D)),B_range[1])

    path = "2DCrystalLattice/"
    foldername = "$(uniqueID)/"
    mkpath(path*foldername)

    #Save params of B field range for plotting
    h5open("$(path)$(foldername)params$(convert(Int64,round(startB))).h5","w") do fid
        current_B = create_group(fid,"params")
        current_B["Bs"] = Vector{Float64}(B_range)
    end

    filepath = "$(path)$(foldername)$(campaign).h5"
    header = ["Current","Electric Field"]
    h5open(filepath,"w") do fid
        #put simulation grid for campaign in file
        sim_grid = create_group(fid,"simulation_grid")
        sim_grid["α"] = start_α
        sim_grid["β"] = start_β
        sim_grid["1_m"] = start_m⁻¹
        sim_grid["σ"] = start_σ
        
        #append metadata for whole campaign to grid
        for (key,val) in metadata 
            HDF5.attributes(sim_grid)[key] = val
        end
        println("Setup Complete")

        #create folder for current campaign
        campaign_group = create_group(fid,"data")
    
        #iterate through B fields, recording data and shot-specific metadata
        for B in B_range
            finder, α, β, m⁻¹,σ = simulation_setup(
            vortex_radius,N,num_crystal,thickness,
            tstep,GL,init_σ,norm_resist,norm_inv_mass,
            Ecrit,Jramp,holdtime,init_hold,xmin,ymin,
            yperiodic,alphaN,betaN,FindType,levelcount,
            tol,backend,string(pkgversion(TDGL2D)),
            B)

            println("Running simulation with B = $(B)")
            sim_data, timetaken = find_jc(finder)

            shot_metadata = Dict("Applied field" => B,"Time taken" => timetaken)
            data_group = create_group(campaign_group,"$(B)field_data")
            for (i,h) in enumerate(header)
                data_group["$h"] = sim_data[i]
            end
            for (key,val) in shot_metadata 
                HDF5.attributes(data_group)[key] = val
            end
        end 
    end 
    println("Simulation complete, time taken = $(time()-init_time)")  
end

function main(uID="default",st=95,inc=10,stp=10,vortex_radius=4,tstep=0.5,GL=10.0,levelcount=3,tol=1e-3,init_σ=1.0,norm_resist=2.0,norm_inv_mass=2.0,Ecrit=10.0^(-4),Jramp=10.0^(-5),holdtime=500,N=2,num_crystal=0,thickness=1.0,xmin=1,ymin=1,yperiodic=true,alphaN=-1,betaN=1,backend="CUDA")
    run_simulation(uID,st,st+inc,stp,vortex_radius,tstep,GL,levelcount,tol,init_σ,norm_resist,norm_inv_mass,Ecrit,Jramp,holdtime,N,num_crystal,thickness,xmin,ymin,yperiodic,alphaN,betaN,backend)
end

if length(ARGS) == 0 
    println("Started with default values")
    main()
else
    main(ARGS[1],
    parse(Float64,ARGS[2]),
    parse(Float64,ARGS[3]),
    parse(Float64,ARGS[4]),
    parse(Int64,ARGS[5]),
    parse(Float64,ARGS[6]),
    parse(Float64,ARGS[7]),
    parse(Int64,ARGS[8]),
    parse(Float64,ARGS[9]),
    parse(Float64,ARGS[10]),
    parse(Float64,ARGS[11]),
    parse(Float64,ARGS[12]),
    parse(Float64,ARGS[13]),
    parse(Float64,ARGS[14]),
    parse(Int64,ARGS[15]),
    parse(Int64,ARGS[16]),
    parse(Int64,ARGS[17]),
    parse(Float64,ARGS[18]),
    parse(Int64,ARGS[19]),
    parse(Int64,ARGS[20]),
    parse(Bool,ARGS[21]),
    parse(Float64,ARGS[22]),
    parse(Float64,ARGS[23]),
    ARGS[24]) 
end


