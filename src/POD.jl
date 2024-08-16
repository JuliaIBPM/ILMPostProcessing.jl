struct PODModes{DT}
    Xmean::DT
    Xnorm::Vector{DT}
    phi::Vector{DT}
    a::Matrix{Float64}
    lambda::Vector{Float64}
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
"""
function PODModes(X::Vector{T}; tolerance=0.99) where T
    Xmean = mean(X)
    Xnorm = map(col -> col - Xmean, X) # normalized by mean

    # calculate X transpose * X matrix and find its eigenvectors/values
    XTX = [dot(xi,xj) for xi in Xnorm, xj in Xnorm]
    lambda, psi = eigen(XTX,sortby=-) # sorts from largest to smallest

    # filter out eigenvalues based on energy tolerance
    lambda_cumsum = cumsum(lambda)
    r = findfirst(lambda_cumsum .>= tolerance*lambda_cumsum[end])
    # perform truncation of modes
    lambda_trunc = lambda[1:r]
    psi_trunc = psi[:,1:r]


    # calculate POD modes
    #phi = [mapreduce((Xi,psi_ij) -> Xi .* psi_ij/sqrt(lambda_i), +, X, psicol) for (psicol,lambda_i) in zip(eachcol(psi_trunc), lambda_trunc)] 
    phi = _calculate_U(X,psi_trunc,sqrt.(lambda_trunc))
    a = [dot(Xk, phi_j) for Xk in Xnorm, phi_j in phi]

    # reconstructed flow field at last solved timestep, ensuring mean is added back
    # fieldReconst = mapreduce((aj, phi_j) -> aj .* phi_j, +, a[end,:], phi) + Xmean
    # return PODModes{typeof(Xnorm[1]),typeof(fieldReconst)}(Xnorm, phi, a, fieldReconst)
    return PODModes{typeof(Xmean)}(Xmean, Xnorm, phi, a, lambda_trunc)
end