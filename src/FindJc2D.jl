const EXPONENTIAL = true
const SUBSECTION = false

"wrapper around the london multigrid solver for determining the critical current density in 2d systems"
mutable struct JC2DFinder{R,VR,VC} <: Finder
    solver::ImplicitLondonMultigridSolver{2,R,RectPrimalForm0Data{2,R,VR},RectPrimalForm1Data{2,R,VR},RectPrimalForm0Data{2,Complex{R},VC},RectPrimalForm1Data{2,Complex{R},VC}}
    mode::JC2DMode
    ecrit::R
    initholdtime::Int64
    jholdtime::Int64
    timesteps::Int64
    curholdsteps::Int64
    j::R
    jrelstep::R
    E_field::R
    esum::R
    δda_rhs::RectPrimalForm1Data{2,R,VR}
    B_field::R
    prev_state::State
    lifetimeE::Vector{R}
end

"constructor for JC2DFinder. Does not initialise BCs"
function JC2DFinder(solver, ecrit::R, initholdtime, jholdtime, jinit::R, jrelstep::R, b_field::R) where {R}
    mode = JC2DInitHold()
    timesteps = 0
    curholdsteps = 0
    e_field = zero(R)
    esum = zero(R)
    δda_rhs = MulTDGL.similar(MulTDGL.state(solver).a)
    data(δda_rhs) .= zero(eltype(data(δda_rhs)))
    lifetimeE = Vector{R}([])

    if EXPONENTIAL
        jinit = jrelstep
        jrelstep = 1.01
    end

    JC2DFinder(solver,
               mode,
               ecrit,
               initholdtime,
               jholdtime,
               timesteps,
               curholdsteps,
               jinit,
               jrelstep,
               e_field,
               esum,
               δda_rhs,
               b_field,
               solver.state_u[1],
               lifetimeE)
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

function increment_J!(finder::JC2DFinder)
    if EXPONENTIAL
        finder.j *= finder.jrelstep
    else
        finder.j += finder.jrelstep
    end
end

function calculate_E!(finder::JC2DFinder,step_data)
    if SUBSECTION
        finder.E_field = E_field_avg(finder.solver,finder.prev_state.a,8)
        finder.prev_state = copy(state(finder))
    else
        finder.E_field = step_data.e[1]
    end
end

"Run simulation for a given B field. Ramp from J_init until E crit is exceeded"
function step!(finder::JC2DFinder)
    sys = system(finder)
    parameters = sys.p

    #record path of E field over time
    push!(finder.lifetimeE,finder.E_field)

    jc2d_bcs!(finder,sys)

    #call london multigrid
    step_data = MulTDGL.step!(finder.solver, finder.δda_rhs, (finder.j,0.0)) 
    calculate_E!(finder,step_data)

    finder.timesteps += 1
    finder.curholdsteps += 1
    finder.esum += finder.E_field

    if finder.mode == JC2DInitHold()
        if finder.curholdsteps > finder.initholdtime / parameters.k - 0.5
            finder.curholdsteps = 0
            finder.mode = JC2DJHold()
            finder.esum = zero(typeof(finder.esum))
        end
    elseif finder.mode == JC2DJHold()
        if finder.jrelstep > zero(typeof(finder.jrelstep))
            if finder.E_field < finder.ecrit
                finder.curholdsteps = 0
                increment_J!(finder)
                finder.esum = zero(typeof(finder.esum))
            elseif finder.curholdsteps > finder.jholdtime / parameters.k - 0.5
                finder.curholdsteps = 0
                finder.mode = JC2DDone()
            end
        else
            if finder.curholdsteps > finder.jholdtime / parameters.k - 0.5
                eavg = finder.esum / finder.curholdsteps
                if eavg < finder.ecrit
                    finder.curholdsteps = 0
                    finder.mode = JC2DDone()
                else
                    finder.curholdsteps = 0
                    increment_J!(finder)   #linear decrease in J ramping
                    finder.esum = zero(typeof(finder.esum))
                end
            end
        end
    end
end

"Run stepping code and record electric field, current and time taken"
function find_jc(f_jc::JC2DFinder,verbose::Bool=true)
    current = Vector{Float64}([])
    b_field = Vector{Float64}([])
    mode = Vector{String}([])

    starttime = time()
    while f_jc.mode != JC2DDone()
        push!(b_field,f_jc.B_field)
        push!(mode,string(f_jc.mode))
        push!(current,f_jc.j)
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

    return [current,f_jc.lifetimeE,b_field,mode], timetaken
end