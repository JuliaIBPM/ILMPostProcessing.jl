module ILMPostProcessing

using LinearAlgebra
using Statistics
using OrdinaryDiffEq
using Interpolations
using RecursiveArrayTools
using HCubature

export pod, dmd, PODModes, DMDModes, compute_FTLE!, compute_LAVD!, compute_IVD!, compute_trajectory,
        compute_streamline, compute_streakline, displacement_field, Trajectories,
        field_along_trajectory, VectorFieldSequence, ScalarFieldSequence



include("utilities.jl")
include("trajectories.jl")
include("POD.jl")
include("DMD.jl")
include("FTLE.jl")
include("LAVD.jl")

include("plot_recipes.jl")


end # module ILMPostProcessing
