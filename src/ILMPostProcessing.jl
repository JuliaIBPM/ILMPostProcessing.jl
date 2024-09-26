module ILMPostProcessing

using LinearAlgebra
using Statistics
using OrdinaryDiffEq
using Interpolations
using RecursiveArrayTools

export pod, dmd, PODModes, DMDModes, compute_FTLE!, compute_trajectory,
        compute_streamline, compute_streakline, displacement_field, Trajectories,
        field_along_trajectory

include("utilities.jl")
include("trajectories.jl")
include("POD.jl")
include("DMD.jl")
include("FTLE.jl")

include("plot_recipes.jl")


end # module ILMPostProcessing
