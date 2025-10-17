module TDGL_Polycrystal

using MulTDGL
using LinearAlgebra
using StaticArrays
using KernelAbstractions
using Dates
using HDF5
using Adapt
using ArgParse
using CUDA
using Random

const Version = "0.1.7" 

include("Defs.jl")
include("Utilities.jl")
include("Voronoi.jl")
include("Crystal2D.jl")
include("FindJc2D.jl")
include("FindBFix.jl")
include("FindLinX.jl")
include("FindEvsJ.jl")
include("Setup.jl")
include("CLI.jl")

export JC2DFinder, Bfixed, BVarLinXFinder, simulation_setup, parse_CL,
       to_string, save_metadata, save_simdata, find_jc, apply_pattern,
       JC2DInitHold, JC2DJHold, JC2DDone, BVarLinX, EvsJ, step!
end #module TDGL_Polycrystal
