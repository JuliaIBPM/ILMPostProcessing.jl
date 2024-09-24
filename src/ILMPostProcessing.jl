module ILMPostProcessing

using ImmersedLayers
using LinearAlgebra
using Statistics

export pod, dmd, PODModes, DMDModes, make_interp_fields!, gen_init_conds, euler_forward, euler_backward, adams_bashforth_2_forward, adams_bashforth_2_backward, compute_FTLE!

include("utilities.jl")
include("POD.jl")
include("DMD.jl")
include("FTLE.jl")


end # module ILMPostProcessing
