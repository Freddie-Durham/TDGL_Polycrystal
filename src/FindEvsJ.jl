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
    δda_rhs::RectPrimalForm1Data{2,R,VR}
end

function EvsJ(solver, ecrit::R, shortholdtime, longholdtime, jinit::R, jrelstep::R, startB::R, max_steps::Int64) where {R}
    mode = JC2DInitHold()
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
    
    #update E field
    finder.E_field = step_data.e[1]
    #calculate supercurrent density
    finder.js = collect(Js_avg(finder.solver, sys, st))

    finder.timesteps += 1
    finder.curholdsteps += 1

    if finder.j > finder.jrelstep * finder.max_steps
        finder.mode = JC2DDone()
    elseif finder.mode == JC2DInitHold() && finder.curholdsteps >= finder.longholdtime / parameters.k
        finder.curholdsteps = 0
        finder.mode = JC2DJHold()
    elseif finder.mode == JC2DJHold() && finder.curholdsteps >= finder.shortholdtime / parameters.k
        finder.j += finder.jrelstep
        finder.curholdsteps = 0
    end
end

function find_jc(f_jc::EvsJ,verbose::Bool=true)
    current = Vector{Float64}([])
    super_current = Vector{Float64}([])
    b_field = Vector{Float64}([])
    e_field = Vector{Float64}([])

    starttime = time()
    while f_jc.mode != JC2DDone()
        push!(b_field,f_jc.B_field)
        push!(current,f_jc.j)
        push!(super_current,f_jc.js[1])
        push!(e_field,f_jc.E_field)
        if verbose && (f_jc.curholdsteps % f_jc.shortholdtime == 0)
            @time step!(f_jc)
            println("Time taken = $(time()-starttime)")
            println("Current = $(f_jc.j)")
            println("Supercurrent = $(f_jc.js)")
            println("Electric Field = $(f_jc.E_field)")
            println("Magnetic Field = $(f_jc.B_field)")
            println("Current hold: $(f_jc.curholdsteps)")
        else
            step!(f_jc)
        end
    end
    timetaken = time()-starttime

    return [current,super_current,b_field,e_field], timetaken
end