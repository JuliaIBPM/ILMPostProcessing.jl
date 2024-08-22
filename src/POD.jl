struct PODModes{DT,AT}
    Xmean::DT
    Xnorm::Vector{DT}
    phi::Vector{DT}
    a::Matrix{AT}
    lambda::Vector{AT}
end

"""
    PODModes(X::Vector{T}[; tolerance=0.99])

Calculate the POD modes and associated time-varying coefficients from an array
of snapshot data `X`. This `X` plays the role of a snapshot matrix, whose columns are
snapshots of the data. However, it is actually to be stored as a type `Vector{T}` where `T<:GridData`.
It can be generated with a function like `velocity(sol,sys,trange)`, where `sol` is a `ODESolution`,
`sys` is an `ILMSystem`, and `trange` is an array of times, e.g., `trange=range(1,10,length=100)`.
The number of POD modes retained in the decomposition is set by `tolerance`: this specifies
the fraction of the total energy to keep, and defaults to 99 percent. 

The output of `PODModes` is a structure with the following fields
- `Xmean`: temporal mean of the data. type `T`
- `Xnorm`: original `X` vector with mean removed. Each element is of type `T`
- `phi`: vector of POD modes. Each element is of type `T`
- `a`: matrix of POD coefficients. Number of columns is same as number of entries in `phi`. Column `k` constitutes the time-varying coefficient for mode `k` in `phi`.
- `lambda`: vector of modal energies, arranged in decreasing order, corresponding to the modes in `phi`
- `psi`: matrix of 
"""
function PODModes(X::AbstractVector{T}; tolerance=0.99) where T
    
    Xmean, Xnorm = _pod_preliminaries(X)

    lambda, psi = _eigen_correlation_matrix(Xnorm)

    # filter out eigenvalues based on energy tolerance
    lambda_cumsum = cumsum(lambda)
    r = findfirst(lambda_cumsum .>= tolerance*lambda_cumsum[end])

    # perform truncation of modes
    lambda_trunc = lambda[1:r]
    psi_trunc = psi[:,1:r]

    _podmodes(Xmean,Xnorm,lambda_trunc,psi_trunc)

end

"""
    PODModes(X::Vector{T},r::Int)

Perform POD on snapshot data `X` and truncate to `r` modes
"""
function PODModes(X::AbstractVector{T},r::Int) where T

    Xmean, Xnorm = _pod_preliminaries(X)

    lambda, psi = _eigen_correlation_matrix(Xnorm)

    # perform truncation of modes
    lambda_trunc = lambda[1:r]
    psi_trunc = psi[:,1:r]

    _podmodes(Xmean,Xnorm,lambda_trunc,psi_trunc)

end

function _pod_preliminaries(X)

    Xmean = mean(X)
    Xnorm = map(col -> col - Xmean, X) # mean removed

    return Xmean, Xnorm
end

function _eigen_correlation_matrix(Xnorm)

    # calculate X transpose * X matrix and find its eigenvectors/values
    XTX = [dot(xi,xj) for xi in Xnorm, xj in Xnorm]
    lambda, psi = eigen(XTX,sortby=-) # sorts from largest to smallest
    
end


function _podmodes(Xmean,Xnorm,lambda_trunc,psi_trunc)

    # calculate POD modes. Note that Ψ = a/sqrt(Λ)
    #phi = [mapreduce((Xi,psi_ij) -> Xi .* psi_ij/sqrt(lambda_i), +, X, psicol) for (psicol,lambda_i) in zip(eachcol(psi_trunc), lambda_trunc)] 
    phi = _calculate_U(Xnorm,psi_trunc,sqrt.(lambda_trunc)) # Φ = X*Ψ/sqrt(Λ)
    #a = [dot(Xk, phi_j) for Xk in Xnorm, phi_j in phi] # Xᵀ*Φ = sqrt(Λ)*Ψ
    a = psi_trunc*Diagonal(sqrt.(lambda_trunc)) # Xᵀ*Φ = sqrt(Λ)*Ψ

    return PODModes{typeof(Xmean),eltype(a)}(Xmean, Xnorm, phi, a, lambda_trunc)

end