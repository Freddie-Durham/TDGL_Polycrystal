module TDGL_Polycrystal

using MulTDGL
using LinearAlgebra
using KernelAbstractions
using Dates
using HDF5
using Adapt
using ArgParse
using CUDA

const Version = "0.1.0" 

include("Defs.jl")
include("Utilities.jl")
include("Crystal2D.jl")
include("Crystal3D.jl")
include("FindJc2D.jl")
include("FindBFix.jl")
include("FindLinX.jl")
include("Setup.jl")
include("CLI.jl")

export JC2D_finder, Bfixed, BVarLinXFinder, simulation_setup, parse_CL
end 