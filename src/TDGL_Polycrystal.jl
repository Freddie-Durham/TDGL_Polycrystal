module TDGL_Polycrystal

using MulTDGL
using LinearAlgebra
using Adapt
using KernelAbstractions
using CUDA
using Dates
using HDF5
using ArgParse

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