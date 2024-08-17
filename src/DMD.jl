struct DMDModes{DT}
    V::Vector{DT}
    mu::Vector{ComplexF64}
end

function DMDModes(Xplus::Vector{T}; tolerance=0.99) where {T}
    nsnap = length(Xplus)
    X = view(Xplus,1:nsnap-1)
    Xp = view(Xplus,2:nsnap)
    pod = PODModes(X;tolerance=tolerance)
    U = pod.phi
    Σ = sqrt.(pod.lambda)
    Ψ = pod.a/Diagonal(Σ)
    Up = _calculate_U(Xp,Ψ,Σ)

    Ã = [dot(ui,upj) for ui in U, upj in Up]
    μ, Ṽ = eigen(Ã)
    
    V = _calculate_U(Up,Ṽ,μ)
    return DMDModes{typeof(V[1])}(V,μ)

end