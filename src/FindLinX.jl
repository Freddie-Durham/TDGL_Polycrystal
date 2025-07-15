"uses Implicit London Multigrid to solve for Jc with a continuous variation in B. Collects a series of E(J) datapoints for extrapolation"
mutable struct BVarLinXFinder{R,VR,VC} <: Finder
    solver::ImplicitLondonMultigridSolver{2,R,RectPrimalForm0Data{2,R,VR},RectPrimalForm1Data{2,R,VR},RectPrimalForm0Data{2,Complex{R},VC},RectPrimalForm1Data{2,Complex{R},VC}}
    mode::LinXMode
    ecrit::R
    shortholdtime::Int64
    longholdtime::Int64
    timesteps::Int64
    curholdsteps::Int64
    j::R
    jrelstep::R
    E_field::R
    B_field::R
    Brelstep::R
    δda_rhs::RectPrimalForm1Data{2,R,VR}
    num_samples::Int64
    sample_index::Int64
    B_index::Int64
    num_Bs::Int64
    allJ::Array{R}
    allE::Array{R}
    Ehist::Vector{R}
    num_holds::Int64
    max_holds::Int64
end

"constructor for BVarLinXFinder. Does not initialise BCs"
function BVarLinXFinder(solver, ecrit::R, shortholdtime, longholdtime, jinit::R, jrelstep::R, startB::R, stopB::R, stepB::R, num_samples) where {R}
    mode = JC2DInitHold()
    timesteps = 0
    curholdsteps = 0
    num_holds = 0
    max_holds = 10 #set how many times we try to wait for true equilibrium
    e_field = zero(R)
    Ehist = Vector{R}([])
    δda_rhs = MulTDGL.similar(MulTDGL.state(solver).a)
    data(δda_rhs) .= zero(eltype(data(δda_rhs)))

    #ensure user input error doesn't lead to solver running indefinitely
    if startB > stopB 
        stepB = -abs(stepB)
    else
        stepB = abs(stepB)
    end
    @assert stepB != zero(R)

    sample_index = 2
    B_index = 1
    num_Bs = convert(Int64,1+round(abs((startB-stopB)/stepB)))
    #1st index of each row contains B field
    allJ = Array{R}(undef,num_Bs,num_samples+1)
    allE = Array{R}(undef,num_Bs,num_samples+1)
    allE[B_index,1] = startB

    BVarLinXFinder(solver,
               mode,
               ecrit,
               shortholdtime,
               longholdtime,
               timesteps,
               curholdsteps,
               jinit,
               abs(jrelstep),
               e_field,
               startB,
               stepB,
               δda_rhs,
               num_samples,
               sample_index,
               B_index,
               num_Bs,
               allJ,
               allE,
               Ehist,
               num_holds,
               max_holds)
end


function reached_equilibrium!(finder::BVarLinXFinder,tol=0.05)
    #record E(J) at equilibrium then increment/extrapolate J or move to next B
    finder.curholdsteps = 0
    Eav = period_avg(finder.Ehist,2*finder.shortholdtime)
    Eupper = period_avg(finder.Ehist,finder.shortholdtime)
    Elower = period_avg(finder.Ehist[1:end-finder.shortholdtime],finder.shortholdtime)
    
    #If avg E contains a significant trend, we continue holding
    if (isapprx(abs(Eupper-Elower)+Eav, Eav, tol/2) & (Eav > 0)) || (finder.num_holds >= finder.max_holds) 


        Eav = abs(Eav)       #prevent crashes if Eav<0 after longholdtime*max_holds sim-steps
        finder.num_holds = 0
        finder.Ehist = Vector{typeof(finder.E_field)}([])
        finder.allJ[finder.B_index,finder.sample_index] = finder.j
        finder.allE[finder.B_index,finder.sample_index] = Eav

        curr_E_vals = finder.allE[finder.B_index,2:finder.sample_index]

        #after we ramp up to Ecrit the first time, we ramp one more time
        #this lets us determine our first estimate for m
        if (finder.sample_index == 2) && (finder.B_index == 1)
            finder.j = finder.j - finder.jrelstep  
            finder.sample_index += 1

        #after collecting num_samples of data, move on to next B field
        #1st index of allE stores B field, 1st index of allJ stores m
        elseif (finder.sample_index >= finder.num_samples + 1) || valfound(curr_E_vals,finder.ecrit,tol)
            #special case if very first two data points are accurate within tolerance, still need to record gradient
            if (finder.sample_index == 3) & (finder.B_index == 1)
                m = lin_ext(log.(finder.allJ[finder.B_index,2:finder.sample_index]),log.(curr_E_vals))
                finder.allJ[finder.B_index,1] = m  
            end
            
            #we dont need the rest of this row as we have found Ec within tolerance
            finder.allE[finder.B_index,finder.sample_index+1:end] .= -1

            finder.sample_index = 2
            finder.B_field += finder.Brelstep #in case of periodic BCs, may need to change c and dc according to new flux
            finder.B_index += 1
            if finder.B_index <= finder.num_Bs
                finder.allE[finder.B_index,1] = finder.B_field
            else
                finder.mode = JC2DDone()
            end

        #When we collect 2 data points, every subsequent data point is found by extrapolation to Ecrit
        #We can guess the 2nd data point using the gradient found at the previous B field 
        else
            m = 0
            if finder.sample_index == 2
                m = finder.allJ[finder.B_index-1,1]
            else
                newm = lin_ext(
                    log.(finder.allJ[finder.B_index,2:finder.sample_index]),
                    log.(curr_E_vals))
                oldm = finder.allJ[finder.B_index,1]
                #take average of current predicted value and previous predicted value of gradient
                m = (newm + oldm) / 2
            end
            #record new gradient for future reference
            finder.allJ[finder.B_index,1] = m   

            #adjust Ecrit to find values slightly above and below, within tolerance
            Etarget = finder.ecrit
            count_above = count([is_above(E,finder.ecrit,tol) for E in curr_E_vals])
            count_below = count([is_below(E,finder.ecrit,tol) for E in curr_E_vals])
            
            Etarget *= (1-tol/2)^count_above
            Etarget *= (1+tol/2)^count_below

            newj = exp(invert_linear(m,log(finder.j),log(Eav),log(Etarget)))
            #prevent values of J that are way too high (hard limits to ensure simulation doesn't break)
            if isnan(newj) | isinf(newj)
                println(newj)
                newj = finder.allJ[finder.B_index,finder.sample_index-1]+finder.jrelstep
            end
            finder.j = max(min(newj,5e-2),1e-7)
            finder.sample_index += 1
        end
    else
        finder.num_holds += 1
    end
end

"Run simulation for a given range of B fields. Collect num_samples of E(J) data"
function step!(finder::BVarLinXFinder)
    sys = system(finder)
    parameters = sys.p

    #record path of E field over time
    push!(finder.Ehist,finder.E_field)

    #update boundary conditions
    jc2d_bcs!(finder,sys)

    #call london multigrid
    step_data = MulTDGL.step!(finder.solver, finder.δda_rhs, (finder.j,0.0)) 
    
    finder.E_field = step_data.e[1]

    finder.timesteps += 1
    finder.curholdsteps += 1

    #If in initial hold state, do nothing unless we have exceeded our stated hold time.
    #Then move to a state of being above or below Jc (this is a rough estimate based on requested Ecrit)
    #We re-enter this state every time we increase the B field
    if finder.mode == JC2DInitHold()
        if finder.curholdsteps >= finder.shortholdtime / parameters.k
            finder.curholdsteps = 0
            finder.mode = JC2DJHold()
        end
    elseif finder.mode == JC2DJHold()
        if finder.E_field < finder.ecrit
            finder.j += finder.jrelstep
            finder.curholdsteps = 0
        elseif finder.curholdsteps >= finder.longholdtime / parameters.k
            reached_equilibrium!(finder)
            finder.mode = BVarLinX()
        end
    #If we are below Jc, we can wait until the E field finds new equilibrium 
    #We record E(J) at equilibrium then step up J_relstep
    #Repeat num_samples times, recording B field and E(J) data in allJ,allE
    #Then increase B and go back into inithold
    elseif finder.mode == BVarLinX()
        if finder.curholdsteps >= finder.longholdtime / parameters.k
            reached_equilibrium!(finder)
        end
    end
end


"Run LinX BVar code and record array of E(J) values for a range of B fields"
function find_jc(f_jc::BVarLinXFinder,verbose::Bool=true)
    current = Vector{Float64}([])
    b_field = Vector{Float64}([])
    e_field = Vector{Float64}([])

    starttime = time()
    while f_jc.mode != JC2DDone()
        push!(b_field,f_jc.B_field)
        push!(current,f_jc.j)
        push!(e_field,f_jc.E_field)
        if verbose
            @time step!(f_jc)
            println("Time taken = $(time()-starttime)")
            println("Current = $(f_jc.j)")
            println("Electric Field = $(f_jc.E_field)")
            println("Magnetic Field = $(f_jc.B_field)")
            println("Mode = "*string(f_jc.mode))
            println("Current hold: $(f_jc.curholdsteps)")
            println("Hold Attempts: $(f_jc.num_holds)")
        else
            step!(f_jc)
        end
    end
    timetaken = time()-starttime

    return [current,b_field,e_field,f_jc.allJ,f_jc.allE], timetaken
end