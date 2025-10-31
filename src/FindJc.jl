const EXPONENTIAL = true
const SUBSECTION = false

"wrapper around the london multigrid solver for determining the critical current density"
mutable struct JcFinder{N,R,VR,VC} <: Finder
    solver::ImplicitLondonMultigridSolver{N,R,RectPrimalForm0Data{N,R,VR},RectPrimalForm1Data{N,R,VR},RectPrimalForm0Data{N,Complex{R},VC},RectPrimalForm1Data{N,Complex{R},VC}}
    mode::JcMode
    ecrit::R
    initholdtime::Int64
    jholdtime::Int64
    timesteps::Int64
    curholdsteps::Int64
    j::R
    jrelstep::R
    E_field::R
    esum::R
    δda_rhs::RectPrimalForm1Data{N,R,VR}
    B_field::R
    prev_state::State
end

"constructor for JcFinder. Does not initialise BCs"
function JcFinder(solver::ImplicitLondonMultigridSolver{N,R,VR,VC}, ecrit, initholdtime, jholdtime, jinit, jrelstep, b_field) where {N,R,VR,VC}
    mode = JcInitHold()
    timesteps = 0
    curholdsteps = 0
    e_field = zero(R)
    esum = zero(R)
    δda_rhs = MulTDGL.similar(MulTDGL.state(solver).a)
    data(δda_rhs) .= zero(eltype(data(δda_rhs)))

    JcFinder(solver,
               mode,
               R(ecrit),
               initholdtime,
               jholdtime,
               timesteps,
               curholdsteps,
               R(jinit),
               R(jrelstep),
               e_field,
               esum,
               δda_rhs,
               R(b_field),
               solver.state_u[1])
end

"find the average value of a horizontal fraction of an array (centered around the midpoint) then divide by the ratio of the fraction to the total area"
function subarray_avg(sarr,extent,frac) 
    y_slice = frac == 0 ? 0 : extent ÷ frac
    avg = sum(asarray(sarr,1)[:,1+y_slice:end-y_slice])
    return avg / (1-2/frac)
end

"calculate the electric field on each edge in the system, then return the average value in the x direction"
function E_field_avg(s::ImplicitLondonMultigridSolver,a_prev,frac=0)
    m = mesh(s)
    δt = parameters(s).k
    δh = m.h[1]

    set_form!((e,_,a,a_prev) -> a[e]-a_prev[e], s.scratch_1, MulTDGL.state(s).a, a_prev) 
    dAdt_avg = subarray_avg(s.scratch_1,m.extent[2],frac) /δt

    set_form!((e,_,ϕ) -> dₑ(m,ϕ,e), s.scratch_1, s.state_u[1].φ)
    ∇ϕ_avg = subarray_avg(s.scratch_1,m.extent[2],frac) /δh

    return (-dAdt_avg -∇ϕ_avg) * δh / measure(m)
end

function increment_J!(finder::JcFinder)
    if EXPONENTIAL
        finder.j *= finder.jrelstep
    else
        finder.j += finder.jrelstep
    end
end

function calculate_E!(finder::JcFinder,step_data)
    if SUBSECTION
        finder.E_field = E_field_avg(finder.solver,finder.prev_state.a,8)
        finder.prev_state = copy(state(finder))
    else
        finder.E_field = step_data.e[1]
    end
end

"Run simulation for a given B field. Ramp from J_init until E crit is exceeded"
function step!(finder::JcFinder{N}) where {N}
    sys = system(finder)
    parameters = sys.p

    jc_bcs!(finder,sys)
    current_array = zeros(eltype(finder.j),N)
    current_array[1] = finder.j
    applied_current = NTuple{N,eltype(finder.j)}(current_array)

    #call london multigrid
    step_data = MulTDGL.step!(finder.solver, finder.δda_rhs, applied_current) 
    
    calculate_E!(finder,step_data)
    finder.timesteps += 1
    finder.curholdsteps += 1
    finder.esum += finder.E_field

    if finder.mode == JcInitHold()
        if finder.curholdsteps > finder.initholdtime / parameters.k - 0.5
            finder.curholdsteps = 0
            finder.mode = JcJHold()
            finder.esum = zero(typeof(finder.esum))
        end
    elseif finder.mode == JcJHold()
        if finder.jrelstep > zero(typeof(finder.jrelstep))
            if finder.E_field < finder.ecrit
                finder.curholdsteps = 0
                increment_J!(finder)
                finder.esum = zero(typeof(finder.esum))
            elseif finder.curholdsteps > finder.jholdtime / parameters.k - 0.5
                finder.curholdsteps = 0
                finder.mode = JcDone()
            end
        else
            if finder.curholdsteps > finder.jholdtime / parameters.k - 0.5
                eavg = finder.esum / finder.curholdsteps
                if eavg < finder.ecrit
                    finder.curholdsteps = 0
                    finder.mode = JcDone()
                else
                    finder.curholdsteps = 0
                    increment_J!(finder)   
                    finder.esum = zero(typeof(finder.esum))
                end
            end
        end
    end
end

"Run stepping code and record electric field, current and time taken"
function find_jc(f_jc::JcFinder,verbose::Bool=true)
    current = Vector{Float64}([])
    b_field = Vector{Float64}([])
    e_field = Vector{Float64}([])
    mode = Vector{String}([])

    starttime = time()
    while f_jc.mode != JcDone()
        push!(b_field,f_jc.B_field)
        push!(mode,string(f_jc.mode))
        push!(current,f_jc.j)
        push!(e_field,f_jc.E_field)
        if verbose
            @time step!(f_jc)
            println("Time taken = $(time()-starttime)")
            println("Current = $(f_jc.j)")
            println("Electric Field = $(f_jc.E_field)")
            println("Magnetic Field = $(f_jc.B_field)")
            println("Mode = "*string(f_jc.mode))
        else
            step!(f_jc)
        end
    end
    timetaken = time()-starttime

    return [current,e_field,b_field,mode], timetaken
end