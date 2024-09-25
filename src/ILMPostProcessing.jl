module ILMPostProcessing

using ImmersedLayers
using LinearAlgebra
using Statistics
using OrdinaryDiffEq

export pod, dmd, PODModes, DMDModes, compute_FTLE!, compute_trajectory,
        compute_streamline, compute_streakline, displacement_field, Trajectories,
        pick_trajectory

include("utilities.jl")
include("trajectories.jl")
include("POD.jl")
include("DMD.jl")
include("FTLE.jl")


end # module ILMPostProcessing
