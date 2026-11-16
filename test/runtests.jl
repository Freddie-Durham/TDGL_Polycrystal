using TDGL_Polycrystal
using MulTDGL
using KernelAbstractions
using Adapt
using CUDA
using HDF5
using Test

const B_app = 0.5
const pixels = 1
const tstep = 0.5
const GL = 10.0
const init_σ = 1.0
const norm_res = 2.0
const norm_mass = 2.0
const Ecrit = 1e-4
const Jramp = 1e-4
const init_hold_time = 2
const wait_time = 2
const xmin = 8
const ymin = 8
const zmin = 8
const alphaN = -1.0
const betaN = 1.0
const levelcount = 2
const tol = 1e-3
const yperiodic = true

@testset "TDGL Polycrystal Tests" begin
    @testset "utilities.jl" begin
        include("test_utils.jl")
    end

    @testset "voronoi.jl" begin
        include("test_voronoi.jl")
    end

    @testset "pattern.jl" begin
        include("test_pattern.jl")
    end

    @testset "setup.jl" begin
        include("test_setup.jl")
    end

    @testset "find_jc.jl" begin
        include("test_findjc.jl")
    end

    @testset "cli.jl" begin
        include("test_cli.jl")
    end
end
