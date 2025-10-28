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

const Version = "0.1.7" 

include("Defs.jl")
include("Utilities.jl")
include("Voronoi.jl")
include("Crystal2D.jl")
include("FindJc.jl")
include("FindBFix.jl")
include("FindLinX.jl")
include("FindEvsJ.jl")
include("Setup.jl")
include("CLI.jl")

export JcFinder, Bfixed, BVarLinXFinder, simulation_setup, simulation_setup_3D,
    parse_CL, to_string, save_metadata, save_simdata, find_jc, apply_pattern,
    JcInitHold, JcJHold, JcDone, BVarLinX, EvsJ, step!
end #module TDGL_Polycrystal
