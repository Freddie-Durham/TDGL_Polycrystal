module TDGL_Polycrystal

using CUDA
using StaticArrays
using MulTDGL
using LinearAlgebra
using KernelAbstractions
using Dates
using HDF5
using Adapt
using ArgParse
using Random

const Version = "0.1.9" 

include("utilities.jl")
include("voronoi.jl")
include("pattern.jl")
include("find_jc.jl")
include("setup.jl")
include("cli.jl")

export JcFinder, simulation_setup, simulation_setup_3D,
    parse_CL, save_metadata, save_simdata, find_jc, apply_pattern,
    JcInitHold, JcJHold, JcDone, step!, apply_B_field!
end #module TDGL_Polycrystal
