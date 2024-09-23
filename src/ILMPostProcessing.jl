module ILMPostProcessing

using ImmersedLayers
using LinearAlgebra
using Statistics

export pod, dmd, PODModes, DMDModes, compute_FTLE!, compute_trajectory

include("utilities.jl")
include("trajectories.jl")
include("POD.jl")
include("DMD.jl")
include("FTLE.jl")


end # module ILMPostProcessing
