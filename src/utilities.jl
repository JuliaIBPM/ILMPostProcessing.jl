"""
    ddt(u::AbstractVector,Δt[,mydiff=:backward_diff])

 Calculate the time derivative of vector data `u`, with time step size `Δt`.
 The default method is backward differencing, but this can be changed to
 `:forward_diff` or `:central_diff`.
 """
 function ddt(u::AbstractVector{T},t::AbstractVector{S};mydiff::Symbol=:forward_diff) where {T,S}
    du = eval(mydiff)(u)./eval(mydiff)(t)
    return du
 end

 # Some basic differencing routines
 function backward_diff(u::AbstractVector{T}) where {T}
     du = zero(u)
     du[2:end] .= u[2:end] .- u[1:end-1]
     du[1] = du[2]
     return du
 end
 function forward_diff(u::AbstractVector{T}) where {T}
     du = zero(u)
     du[1:end-1] .= u[2:end] .- u[1:end-1]
     du[end] = du[end-1]
     return du
 end
 function central_diff(u::AbstractVector{T}) where {T}
     du = zero(u)
     du[2:end-1] .= 0.5*u[3:end] .- 0.5*u[1:end-2]
     du[1] = du[2]
     du[end] = du[end-1]
     return du
 end

