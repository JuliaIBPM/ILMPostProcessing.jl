using Documenter, ILMPostProcessing

ENV["GKSwstype"] = "nul" # removes GKS warnings during plotting

makedocs(
    sitename = "ILMPostProcessing.jl",
    doctest = true,
    clean = true,
    modules = [ILMPostProcessing],
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Manual" => ["manual/pod.md",
                     "manual/dmdtest.md",
                     "manual/ftle.md",
                     "manual/functions.md"
                     ]
        #"Internals" => [ "internals/properties.md"]
    ],
    #format = Documenter.HTML(assets = ["assets/custom.css"])
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine = MathJax(Dict(
            :TeX => Dict(
                :equationNumbers => Dict(:autoNumber => "AMS"),
                :Macros => Dict()
            )
        ))
    ),
    #assets = ["assets/custom.css"],
    #strict = true
)


#if "DOCUMENTER_KEY" in keys(ENV)
deploydocs(
     repo = "github.com/JuliaIBPM/ILMPostProcessing.jl.git",
     target = "build",
     deps = nothing,
     make = nothing
     #versions = "v^"
)
#end
