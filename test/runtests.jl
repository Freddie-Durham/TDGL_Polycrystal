using TDGL_Polycrystal
using Test

const B_app = 0.95
const B_step = 0.01
const Bmax = 0.03
const Bstart = 0.0
const num_samples = 5
const pixels = 2
const N = 2
const grain_size = 30.0
const crystalangle = 30*π/180
const tstep = 1.0
const GL = 10.0
const init_σ = 1.0
const norm_res = 2.0
const norm_inv_mass = 2.0
const Ecrit = 1e-5
const Jramp = 1e-4
const init_hold_time = 50
const wait_time = 50
const xmin = 1
const ymin = 2
const alphaN = -1
const levelcount = 4
const tol = 1e-3
const yperiodic = false
const Version = "0.0.0"

const Verbose = false
const expensive = false

bknd = "CPU"
backend = CPU()
if CUDA.functional()
    bknd = "CUDA"
    backend = CUDABackend()
end

@testset "Utils" begin
    include("test_utils.jl")
end

@testset "2D Crystal" begin
    include("test_2D.jl")
end