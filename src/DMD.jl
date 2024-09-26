struct DMDModes{DT}
    modes::Vector{DT}
    evals::Vector{ComplexF64}
end

"""
    dmd(Xfull,r)

Compute the first `r` DMD modes from the extended snapshot data in `Xfull`.
Both the original and shifted data are drawn from `Xfull`, that is:
`X = Xfull[1:m-1]`
and
`Xp = Xfull[2:m]`

This returns the DMD modes and DMD eigenvalues in a `DMDModes`
structure, with fields `modes` and `evals`.
"""
function dmd(Xplus::Vector{T}, r::Integer) where {T}
    nsnap = length(Xplus)
    X = view(Xplus,1:nsnap-1)
    Xp = view(Xplus,2:nsnap)

    # Calculate eigenvalues and eigenvectors of the total [X; Xp] matrix
    lambda, psi = _eigen_correlation_matrix(X, Xp)
    Q = psi[:,1:r]

    o = ones(r)
    Xhat = _calculate_U(X,Q,o)
    Xphat = _calculate_U(Xp,Q,o)

    lambdaX, VX = _eigen_correlation_matrix(Xhat)
    ΣX = sqrt.(lambdaX)
    UX = _calculate_U(Xhat,VX,ΣX)
    UXp = _calculate_U(Xphat,VX,ΣX)

    Ã = _calculate_XTY_via_dot(UX,UXp)
    μ, Ṽ = eigen(Ã) #,sortby = x -> abs(angle(x)))
    o = ones(length(μ))
    V = _calculate_U(UX,Ṽ,o)

    return DMDModes{typeof(V[1])}(V,μ)

end