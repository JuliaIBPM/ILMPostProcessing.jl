# TRAJECTORY CALCULATION #

const DEFAULT_DT = 0.001
const DEFAULT_DT_STREAK = 0.01

## APIs ##


"""
   compute_trajectory(vel::Edges,sys,X₀::Vector/Vector{Vector},Tmax[,Δt=0.001])

Calculate the trajectory of a particle with initial location `X₀`. The argument
`vel` is edge-type grid data, `sys` is a Navier-Stokes type system, `Tmax` is the final
integration time, and `Δt` is the time step size. The output is the solution
structure for the `OrdinaryDiffEq` package.
"""
compute_trajectory(u::Edges, sys::ILMSystem, X₀::Union{Vector{S},Vector{Vector{S}}},Trange::Tuple;Δt=DEFAULT_DT) where S <: Real =
    compute_trajectory(interpolatable_field(u,sys.base_cache.g)...,X₀,Trange;Δt=Δt)



"""
   compute_trajectory(u,v,X₀::Vector,Trange::Tuple[,Δt=0.001])

Calculate the trajectory of a tracer particle with initial location(s) `X₀`, which
can be specified as either a single vector `[x0,y0]` or a vector of vectors
for multiple tracer particles. The arguments
`u` and `v` are either interpolated velocity field components from a computational solution
or are functions. If they are functions, then each of them should be of the form `u(x,y,t)`
and `v(x,y,t)`; `Trange` is a tuple of the initial and final time of integration; and `Δt` is the
time step size, which defaults to 0.001. The output is the solution
structure for the `OrdinaryDiffEq` package (or, for multiple particles, a vector
of such solution structures).
"""
function compute_trajectory(ufield::AbstractInterpolation{T,2},
                            vfield::AbstractInterpolation{T,2},
                            X₀::Vector{S},Trange::Tuple;Δt::Real=DEFAULT_DT) where {T,S<:Real}

  vfcn!(dR,R,p,t) = _vfcn_autonomous!(dR,R,p,t,ufield,vfield)

  sol = _solve_trajectory(vfcn!,X₀,Trange,Δt)
  return sol

end

function compute_trajectory(ufield::T,vfield::T,
   pts::Vector{Vector{S}},Trange::Tuple;Δt=DEFAULT_DT) where {S<:Real,T<:Union{AbstractInterpolation,Function}}

  sol_array = ODESolution[]
  for X₀ in pts
    sol = compute_trajectory(ufield,vfield,X₀,Trange,Δt=Δt)
    push!(sol_array,sol)
  end
  return sol_array

end

function compute_trajectory(ufcn::Function,
                            vfcn::Function,
                            X₀::Vector{S},Trange::Tuple;Δt::Real=DEFAULT_DT) where {S<:Real}

  velfcn(R,p,t) = _vfcn_nonautonomous(R,p,t,ufcn,vfcn)

  sol = _solve_trajectory(velfcn,X₀,Trange,Δt)
  return sol

end

"""
   compute_streamline(u,v,X₀::Vector,srange::Tuple,t::Real,[,Δt=0.001])

Calculate the streamline(s) passing through location(s) `X₀`, which
can be specified as either a single vector `[x0,y0]` or a vector of vectors
for a rake of streamlines. The arguments
`u` and `v` are either interpolated velocity field components from a computational solution
or are functions. If they are functions, then each of them should be of the form `u(x,y,t)`
and `v(x,y,t)`; `srange` is a tuple of the initial and final time of integration; `t` is
the current time at which the streamline is depicted; and `Δs` is the
time-like step size, which defaults to 0.001. The output is the solution
structure for the `OrdinaryDiffEq` package (or, for multiple points, a vector
of such solution structures).
"""
function compute_streamline(ufcn::Function,
                            vfcn::Function,
                            X₀::Vector{S},srange::Tuple,t::Real;Δs::Real=DEFAULT_DT) where {S<:Real}

  velfcn(R,p,s) = _vfcn_nonautonomous_frozentime(R,p,s,ufcn,vfcn)

  sol = _solve_streamline(velfcn,X₀,srange,Δs,t)
  return sol

end

function compute_streamline(ufield,vfield,
   pts::Vector{Vector{S}},srange::Tuple,t::Real;Δs=DEFAULT_DT) where {S<:Real}

  sol_array = ODESolution[]
  for X₀ in pts
    sol = compute_streamline(ufield,vfield,X₀,srange,t,Δs=Δs)
    push!(sol_array,sol)
  end
  return sol_array

end

"""
   compute_streakline(u,v,X₀::Vector,t[;τmin = t-3.0, Δtstreak=0.01,Δttraj=0.001]) -> Vector, Vector

Calculate a streakline at time `t` for a velocity field `u` and `v`, based on an injection
point `X₀`. The end of the streakline is set by `τmin`, the earliest time
at which a particle passed through the injection point. It defaults to 3 time
units before the current instant `t`. The time step size `Δt` sets the resolution
of the streakline (i.e., how often the particles are sampled along the streakline).
It returns arrays of the x and y coordinates of the streakline.
"""
function compute_streakline(u,v,X₀::Vector{S},t;τmin = t-3.0, Δtstreak::Real=DEFAULT_DT_STREAK, Δttraj::Real=DEFAULT_DT) where {S<:Real}
  τstreak = τmin:Δtstreak:t
  xstreak = zeros(length(τstreak))
  ystreak = zeros(length(τstreak))

  for (i,τ) in enumerate(τstreak)
    traj = compute_trajectory(u,v,X₀,(τ,t),Δt = Δttraj)
    xstreak[i], ystreak[i] = traj.u[end]
  end
  return xstreak, ystreak
end

"""
  field_along_trajectory(f::GridData,sys::NavierStokes,traj::ODESolution[,deriv=0])

Evaluate field `f` (given as grid data) along the trajectory specified by `traj`.
The output is the history of `f` along this trajectory. If `f` is a vector field,
then the component histories are output as a tuple. If `deriv=1`, then it
computes the time derivative of the field along the trajectory. The default
is `deriv=0` (no derivative).
"""
field_along_trajectory(d::GridData,sys,traj;deriv=0) = _field_along_trajectory(d,sys,traj,Val(deriv))




## Internal helper functions ##



function _solve_trajectory(vfcn,u0,Trange,Δt)
  Path = ODEProblem(vfcn,u0,Trange)
  sol = solve(Path,Tsit5(), dt = Δt, maxiters = 1e8, adaptive = false, dense = false)
end

function _solve_streamline(vfcn,u0,Trange,Δt,p)
  Path = ODEProblem(vfcn,u0,Trange,p)
  sol = solve(Path,Tsit5(), dt = Δt, maxiters = 1e8, adaptive = false, dense = false)
end


function _vfcn_autonomous!(dR,R,p,t,u,v)
 dR[1] = u(R[1],R[2])
 dR[2] = v(R[1],R[2])

 return dR
end

function _vfcn_nonautonomous(R,p,t,u,v)
 dR = similar(R)
 dR[1] = u(R[1],R[2],t)
 dR[2] = v(R[1],R[2],t)

 return dR
end

function _vfcn_nonautonomous_frozentime(R,p,t,u,v)
  dR = similar(R)
  dR[1] = u(R[1],R[2],p)
  dR[2] = v(R[1],R[2],p)

 return dR
end

function _field_along_trajectory(v::VectorGridData,sys::ILMSystem,traj::ODESolution,::Val{0})
    vfield_x, vfield_y = interpolatable_field(v,sys.base_cache.g)
   
    vx_traj = eltype(v)[]
    vy_traj = eltype(v)[]
    for x in traj.u
      push!(vx_traj,vfield_x(x...))
      push!(vy_traj,vfield_y(x...))
    end
   
    return vx_traj, vy_traj
end
   
function _field_along_trajectory(s::ScalarGridData,sys::ILMSystem,traj::ODESolution,::Val{0})
    sfield = interpolatable_field(s,sys.base_cache.g)
   
    s_traj = eltype(sfield)[]
    for x in traj.u
      push!(s_traj,sfield(x...))
    end
   
    return s_traj
end
   
function _field_along_trajectory(v::VectorGridData,sys::ILMSystem,traj::ODESolution,::Val{1})
       utraj, vtraj = _field_along_trajectory(v,sys,traj,Val(0))
       return ddt(utraj,traj.t), ddt(vtraj,traj.t)
end
   
_field_along_trajectory(s::ScalarGridData,sys::ILMSystem,traj::ODESolution,::Val{1}) =
       ddt(_field_along_trajectory(s,sys,traj,Val(0)),traj.t)
   