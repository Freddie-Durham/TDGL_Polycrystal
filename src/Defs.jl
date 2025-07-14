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