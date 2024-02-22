struct PODModes{DT, PT}
    Xnorm::Vector{DT}
    phi::Vector{DT}
    a::Matrix{Float64}
    fieldReconst::PT
end

function createSnapshotData(dataFunction::Function, sol, sys::ILMSystem; timestep::Integer=10)
    # extract velocity field at every 10th timestep and normallize by temporal mean
    X = [dataFunction(sol, sys, t) for (u, t) in zip(sol.u[1:timestep:end], sol.t[1:timestep:end])]
end

function PODModes(X::Vector{T}; tolerance=0.99) where T
    Xmean = mean(X)
    Xnorm = map(col -> col - Xmean, X) # normalized by mean

    # calculate X transpose * X matrix and find its eigenvectors/values
    XTX = [dot(xi,xj) for xi in Xnorm, xj in Xnorm]
    F = eigen(XTX)
    lambda = F.values
    psi = F.vectors

    # filter out eigenvalues based on energy tolerance
    lambda_sum = sum(lambda)
    rev_lambda = reverse(lambda)
    lambda_cumsum = cumsum(rev_lambda)
    r = findfirst(lambda_cumsum .>= tolerance*lambda_sum)
    # perform truncation of modes
    lambda_trunc = lambda[length(lambda)-(r-1):end]
    psi_trunc = psi[:,length(lambda)-(r-1):end]

    # calculate POD modes
    phi = [mapreduce((Xi,psi_ij) -> Xi .* psi_ij/sqrt(lambda_i), +, X, psicol) for (psicol,lambda_i) in zip(eachcol(psi_trunc), lambda_trunc)] 
    a = [dot(Xk, phi_j) for Xk in Xnorm, phi_j in phi]

    # reconstructed flow field at last solved timestep, ensuring mean is added back
    fieldReconst = mapreduce((aj, phi_j) -> aj .* phi_j, +, a[end,:], phi) + Xmean
    return PODModes{typeof(Xnorm[1]),typeof(fieldReconst)}(Xnorm, phi, a, fieldReconst)
end