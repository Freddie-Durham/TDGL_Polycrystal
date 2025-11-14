const SUBSECTION = false

abstract type RampMode end

struct Exp_Increase <: RampMode end
struct Exp_Decrease <: RampMode end
struct Linear_Increase <: RampMode end

function decode_method(method::String)
    if uppercase(method) == "EXP_INCREASE"
        return Exp_Increase
    elseif uppercase(method) == "EXP_DECREASE"
        return Exp_Decrease
    elseif uppercase(method) == "LIN_INCREASE"
        return Linear_Increase
    else
        error("Unknown ramp method: $method")
    end
end

struct Ramp_Method{T<:RampMode,R}
    J_init::R
    J_relstep::R
    E_crit::R
    Init_Hold_Time::Int64
    J_Hold_Time::Int64
end

"wrapper around the london multigrid solver for determining the critical current density"
mutable struct JcFinder{N,R,VR,VC,T} <: Finder
    solver::ImplicitLondonMultigridSolver{N,R,RectPrimalForm0Data{N,R,VR},RectPrimalForm1Data{N,R,VR},RectPrimalForm0Data{N,Complex{R},VC},RectPrimalForm1Data{N,Complex{R},VC}}
    mode::JcMode
    curholdsteps::Int64
    ramp_method::Ramp_Method{T}
    j::R
    E_field::R
    esum::R
    δda_rhs::RectPrimalForm1Data{N,R,VR}
    B_field::R
    prev_state::State
end

"constructor for JcFinder. Does not initialise BCs"
function JcFinder(solver::ImplicitLondonMultigridSolver{N,R,VR,VC}, ecrit, initholdtime, jholdtime, jinit, jrelstep, b_field, method) where {N,R,VR,VC}
    mode = JcInitHold()
    curholdsteps = 0
    e_field = zero(R)
    esum = zero(R)
    δda_rhs = MulTDGL.similar(MulTDGL.state(solver).a)
    data(δda_rhs) .= zero(eltype(data(δda_rhs)))

    ramp_method = Ramp_Method{decode_method(method),R}(R(jinit), R(jrelstep), R(ecrit), initholdtime, jholdtime)

    JcFinder(solver,
               mode,
               curholdsteps,
               ramp_method,
               ramp_method.J_init,
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

function next_step!(finder,ramp::Ramp_Method{Exp_Increase,R},parameters) where {R}
    if finder.E_field < ramp.E_crit
        finder.curholdsteps = 0
        finder.j *= R(1 + ramp.J_relstep)
    elseif finder.curholdsteps >= ramp.J_Hold_Time / parameters.k
        finder.curholdsteps = 0
        finder.mode = JcDone()
    end
end

function next_step!(finder,ramp::Ramp_Method{Linear_Increase,R},parameters) where {R}
    if finder.E_field < ramp.E_crit
        finder.curholdsteps = 0
        finder.j += ramp.J_relstep
    elseif finder.curholdsteps >= ramp.J_Hold_Time / parameters.k
        finder.curholdsteps = 0
        finder.mode = JcDone()
    end
end

function next_step!(finder,ramp::Ramp_Method{Exp_Decrease,R},parameters) where {R}
    half_time = round(Int,ramp.J_Hold_Time / (2 * parameters.k))
    
    #start recording E field after half the hold time
    if finder.curholdsteps > half_time
        finder.esum += finder.E_field

        #after full hold time, check average E field vs E crit
        if finder.curholdsteps >= ramp.J_Hold_Time / parameters.k
            if abs(finder.esum) > ramp.E_crit * half_time
                finder.curholdsteps = 0
                finder.j *= R(1 - ramp.J_relstep)
            else
                finder.curholdsteps = 0
                finder.mode = JcDone()
            end
            finder.esum = zero(R)
        end
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
function step!(finder::JcFinder{N,R}) where {N,R}
    sys = system(finder)
    parameters = sys.p

    jc_bcs!(finder,sys)
    current_array = zeros(R,N)
    current_array[1] = finder.j
    applied_current = NTuple{N,R}(current_array)

    #call london multigrid
    step_data = MulTDGL.step!(finder.solver, finder.δda_rhs, applied_current) 
    
    calculate_E!(finder,step_data)
    finder.curholdsteps += 1

    if finder.mode == JcInitHold()
        if finder.curholdsteps >= finder.ramp_method.Init_Hold_Time / parameters.k
            finder.curholdsteps = 0
            finder.mode = JcJHold()
        end
    elseif finder.mode == JcJHold()
       next_step!(finder,finder.ramp_method,parameters)
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