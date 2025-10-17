const Vec3 = SVector{3,Float64}

abstract type Finder end

struct JC2DInitHold end
struct JC2DJHold end
struct JC2DDone end
struct BVarLinX end

Base.string(::JC2DInitHold) = "Initial Hold"
Base.string(::JC2DJHold) = "Hold at fixed J"
Base.string(::JC2DDone) = "Simulation Complete"
Base.string(::BVarLinX) = "Linear fit to find E(J) = Ec"

state(finder::Finder) = MulTDGL.state(finder.solver)
system(finder::Finder) = MulTDGL.system(finder.solver)

const JC2DMode = Union{JC2DInitHold,JC2DJHold,JC2DDone}
const LinXMode = Union{JC2DInitHold,JC2DJHold,BVarLinX,JC2DDone}

"Calculate supercurrent density then return average value"
function Js_avg(solver,sys,st)
    set_form!(solver.scratch_1, sys.m, sys.mat.m⁻¹, sys.u, st) do e, _, m, m⁻¹, u, st
        jₛₑ(m, m⁻¹, u, st.ψ, e)
    end
    MulTDGL.rect_average_in_place!(solver.scratch_1, sys.m)
end