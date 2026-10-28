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

const Version = "0.1.8" 

include("Utilities.jl")
include("Voronoi.jl")
include("Crystal2D.jl")
include("FindJc.jl")
include("Setup.jl")
include("CLI.jl")

export JcFinder, simulation_setup, simulation_setup_3D,
    parse_CL, save_metadata, save_simdata, find_jc, apply_pattern,
    JcInitHold, JcJHold, JcDone, step!, apply_B_field!
end #module TDGL_Polycrystal
