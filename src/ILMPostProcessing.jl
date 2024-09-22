module ILMPostProcessing

using ImmersedLayers
using LinearAlgebra
using Statistics

export pod, dmd, PODModes, DMDModes, compute_FTLE!

include("utilities.jl")
include("POD.jl")
include("DMD.jl")
include("FTLE.jl")


end # module ILMPostProcessing
