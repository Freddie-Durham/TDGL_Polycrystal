module TDGL_Polycrystal

using MulTDGL
using LinearAlgebra
using KernelAbstractions
using Dates
using HDF5
using Adapt
using ArgParse
using CUDA

const Version = "0.1.1" 

include("Defs.jl")
include("Utilities.jl")
include("Crystal2D.jl")
include("Crystal3D.jl")
include("FindJc2D.jl")
include("FindBFix.jl")
include("FindLinX.jl")
include("Setup.jl")
include("CLI.jl")

export JC2DFinder, Bfixed, BVarLinXFinder, simulation_setup, parse_CL,
       to_string, save_metadata, save_simdata, find_jc, tesselate!, tesselateOct!,
       JC2DInitHold, JC2DJHold, JC2DDone, BVarLinX
end 