# ILMPostProcessing.jl

*Tools for post-processing solutions of immersed layer PDEs*

The objective of this package is to supply a variety of post-processing tools for
solutions of PDEs carried out with the [ImmersedLayers.jl](https://github.com/JuliaIBPM/ImmersedLayers.jl) package, and the domain-specific subpackages, such as [ViscousFlow.jl](https://github.com/JuliaIBPM/ViscousFlow.jl). The post-processing tools[^1] currently available
are
* Proper orthogonal decomposition (POD)
* Dynamic mode decomposition (DMD)
* Finite-time Lyapunov exponent (FTLE)
* Lagrangian-averaged vorticity deviation (LAVD)

## Installation

This package works on Julia `1.7` and above and is registered in the general Julia registry. To install from the REPL, type
e.g.,
```julia
] add ILMPostProcessing
```

Then, in any version, type
```julia
julia> using ILMPostProcessing
```

The plots in this documentation are generated using [Plots.jl](http://docs.juliaplots.org/latest/).
You might want to install that, too, to follow the examples.

## References

[^1]: Taira, K. et al (2017) "Modal Analysis of Fluid Flows: An Overview," *AIAA Journal*, 55(12), 4013--4041.
