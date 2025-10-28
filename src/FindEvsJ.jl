const ε0 = 8.854187817e-12 # Vacuum permittivity
const μ0 = 1.256637061e-6 # Vacuum permeability
Jdisp_ratio(ρ,λ) = (ρ^2 * ε0) / (λ^2 * μ0)
const ρ_Nb = 1.5e-7
const λ_Nb = 4.7e-8

const ρ_ReBCO = 2e-6
const λ_ReBCO = 1.5e-7

"Run simulation at a fixed J (≈Jc) and B field to study fluctuations in E field"
mutable struct EvsJ{R,VR,VC} <: Finder
    solver::ImplicitLondonMultigridSolver{2,R,RectPrimalForm0Data{2,R,VR},RectPrimalForm1Data{2,R,VR},RectPrimalForm0Data{2,Complex{R},VC},RectPrimalForm1Data{2,Complex{R},VC}}
    mode::LinXMode
    shortholdtime::Int64  # time for hold so E equlibrates
    longholdtime::Int64   # time for initial hold
    timesteps::Int64
    curholdsteps::Int64 
    max_steps::Int64      # maximum number of Jrel-steps
    j::R
    jrelstep::R
    js::Vector{R}
    E_field::R
    B_field::R
    E_prev::R
    dEdt::R
    Jdisp0::R
    δda_rhs::RectPrimalForm1Data{2,R,VR}
end

function EvsJ(solver, ecrit::R, shortholdtime, longholdtime, jinit::R, jrelstep::R, startB::R, max_steps::Int64) where {R}
    mode = JcInitHold()
    timesteps = 0
    curholdsteps = 0
    δda_rhs = MulTDGL.similar(MulTDGL.state(solver).a)
    data(δda_rhs) .= zero(eltype(data(δda_rhs)))

    EvsJ(solver,
        mode,
        shortholdtime,
        longholdtime,
        timesteps,
        curholdsteps,
        max_steps,
        jinit,
        abs(jrelstep),
        [zero(R),zero(R)],
        zero(R),
        startB,
        zero(R),
        zero(R),
        Jdisp_ratio(ρ_ReBCO,λ_ReBCO),
        δda_rhs,
        )
end

function step!(finder::EvsJ)
    sys = system(finder)
    st = state(finder)
    parameters = sys.p

    #update boundary conditions
    jc2d_bcs!(finder,sys)

    #call london multigrid
    step_data = MulTDGL.step!(finder.solver, finder.δda_rhs, (finder.j,0.0)) 
    
    #update E field and dE/dt
    finder.E_field = step_data.e[1]
    finder.dEdt = (finder.E_field - finder.E_prev) / parameters.k
    finder.E_prev = finder.E_field

    #calculate supercurrent density
    finder.js = collect(Js_avg(finder.solver, sys, st))

    finder.timesteps += 1
    finder.curholdsteps += 1

    if finder.j > finder.jrelstep * finder.max_steps
        finder.mode = JcDone()
    elseif finder.mode == JcInitHold() && finder.curholdsteps >= finder.longholdtime / parameters.k
        finder.curholdsteps = 0
        finder.mode = JcJHold()
    elseif finder.mode == JcJHold() && finder.curholdsteps >= finder.shortholdtime / parameters.k
        finder.j += finder.jrelstep
        finder.curholdsteps = 0
    end
end

function find_jc(f_jc::EvsJ,verbose::Bool=true)
    current = Vector{Float64}([])
    super_current = Vector{Float64}([])
    displace_current = Vector{Float64}([])
    b_field = Vector{Float64}([])
    e_field = Vector{Float64}([])

    starttime = time()
    while f_jc.mode != JcDone()
        J_disp = f_jc.dEdt * f_jc.Jdisp0

        push!(displace_current,J_disp)
        push!(b_field,f_jc.B_field)
        push!(current,f_jc.j)
        push!(super_current,f_jc.js[1])
        push!(e_field,f_jc.E_field)
        if verbose && (f_jc.curholdsteps % f_jc.shortholdtime == 0)
            @time step!(f_jc)
            println("Time taken = $(time()-starttime)")
            println("Current = $(f_jc.j)")
            println("Supercurrent = $(f_jc.js)")
            println("Displacement Current = $(J_disp)")
            println("Electric Field = $(f_jc.E_field)")
            println("Magnetic Field = $(f_jc.B_field)")
            println("Current hold: $(f_jc.curholdsteps)")
        else
            step!(f_jc)
        end
    end
    timetaken = time()-starttime

    return [current,super_current,displace_current,b_field,e_field], timetaken
end