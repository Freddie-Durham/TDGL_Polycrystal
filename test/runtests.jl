using TDGL_Polycrystal
using MulTDGL
using KernelAbstractions
using Adapt
using CUDA
using Test

const B_app = 0.95
const B_step = 0.01
const Bmax = 0.03
const Bstart = 0.0
const num_samples = 5
const pixels = 2
const N = 2
const num_crystal = 1
const grain_thick = 2
const grain_size = 30.0
const crystalangle = 30*π/180
const tstep = 1.0
const GL = 10.0
const init_σ = 1.0
const norm_res = 2.0
const norm_mass = 0.5
const Ecrit = 1e-5
const Jramp = 1e-4
const init_hold_time = 50
const wait_time = 50
const xmin = 32
const ymin = 32
const alphaN = -1
const betaN = 1.0
const levelcount = 4
const tol = 1e-3
const yperiodic = true

const Verbose = false
const expensive = false

bknd = "CPU"
backend = CPU()
if CUDA.functional()
    bknd = "CUDA"
    backend = CUDABackend()
end

@testset "TDGL Polycrystal Tests" begin
    @testset "Utils" begin
        include("test_utils.jl")
    end

    @testset "Crystal Lattice" begin
        include("test_crystal.jl")
    end

    if expensive
        @testset "2D Simulations" begin
            include("test_2D.jl")
        end
    end
end