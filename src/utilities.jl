_calculate_U(X::Vector{T},V::Array,Σ::Vector) where {T} = 
        [mapreduce((Xi,V_ij) -> Xi .* V_ij/σ_i, +, X, Vcol) for (Vcol,σ_i) in zip(eachcol(V), Σ)]