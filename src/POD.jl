struct PODModes{DT}
    Xnorm::Vector{DT}
    phi::Vector{DT}
    a::Matrix{Float64}
end

function createSnapshotData(dataFunction::Function, sol, sys::ILMSystem; timestep::Integer=10)
    # extract velocity field at every 10th timestep and normallize by temporal mean
    X = [dataFunction(sol, sys, t) for (u, t) in zip(sol.u[1:timestep:end], sol.t[1:timestep:end])]
end

function PODModes(X::AbstractVector; tolerance=0.99)
    Xmean = mean(X)
    Xnorm = map(col -> col - Xmean, X) # normalized by mean

    # calculate X transpose * X matrix and find its eigenvectors/values
    XTX = [dot(xi,xj) for xi in Xnorm, xj in Xnorm]
    F = eigen(XTX)
    lambda = F.values
    psi = F.vectors

    # filter out eigenvalues based on energy tolerance
    lambda_sum = cumsum(lambda)
    r = findfirst(lambda_sum .>= tolerance*lambda_sum)

    # calculate POD modes
    phi = [mapreduce((Xi,psi_ij) -> Xi .* psi_ij/sqrt(lambda_i), +, X, psicol) for (psicol,lambda_i) in zip(eachcol(psi), lambda)] 
    a = [dot(Xk, phi_j) for Xk in Xnorm, phi_j in phi]
    return PODModes(Xnorm, phi, a)
end