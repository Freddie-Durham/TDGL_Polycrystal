#alias for 3D Static Vectors
const vec3 = SVector{3,Float64}
#allows broadcast dict like function
(d::Dict)(k) = d[k]

abstract type Finder end

struct JcInitHold end
struct JcJHold end
struct JcDone end

Base.string(::JcInitHold) = "Initial Hold"
Base.string(::JcJHold) = "Hold at fixed J"
Base.string(::JcDone) = "Simulation Complete"

state(finder::Finder) = MulTDGL.state(finder.solver)
system(finder::Finder) = MulTDGL.system(finder.solver)

const JcMode = Union{JcInitHold, JcJHold, JcDone}

"Linear interpolation between min_val and max_val with val between 0 and 1"
linear_interp(min_val, max_val, val) = min_val + (max_val - min_val) * val

MulTDGL_backend(d::MulTDGL.AbstractForm) = KernelAbstractions.get_backend(data(d))

const VERBOSE_LOCK = ReentrantLock()

function verbose_prefix(context)
    task_label = Base.Threads.nthreads() == 1 ? "" : "thread $(Base.Threads.threadid())"

    if isempty(context)
        return isempty(task_label) ? "" : "[$task_label] "
    elseif isempty(task_label)
        return "[$context] "
    else
        return "[$context, $task_label] "
    end
end

"print a complete verbose update without interleaving concurrent threaded output"
function verbose_update(lines; context="")
    lock(VERBOSE_LOCK)
    try
        prefix = verbose_prefix(context)
        for line in lines
            println(prefix, line)
        end
        flush(stdout)
    finally
        unlock(VERBOSE_LOCK)
    end
end
