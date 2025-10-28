#alias for 3D Static Vectors
const vec3 = SVector{3,Float64}
#allows broadcast dict like function
(d::Dict)(k) = d[k]

abstract type Finder end

struct JcInitHold end
struct JcJHold end
struct JcDone end
struct BVarLinX end

Base.string(::JcInitHold) = "Initial Hold"
Base.string(::JcJHold) = "Hold at fixed J"
Base.string(::JcDone) = "Simulation Complete"
Base.string(::BVarLinX) = "Linear fit to find E(J) = Ec"

state(finder::Finder) = MulTDGL.state(finder.solver)
system(finder::Finder) = MulTDGL.system(finder.solver)

const JcMode = Union{JcInitHold,JcJHold,JcDone}
const LinXMode = Union{JcInitHold,JcJHold,BVarLinX,JcDone}

"Calculate supercurrent density then return average value"
function Js_avg(solver,sys,st)
    set_form!(solver.scratch_1, sys.m, sys.mat.m⁻¹, sys.u, st) do e, _, m, m⁻¹, u, st
        jₛₑ(m, m⁻¹, u, st.ψ, e)
    end
    MulTDGL.rect_average_in_place!(solver.scratch_1, sys.m)
end