# TRAJECTORY CALCULATION #

import ImmersedLayers.CartesianGrids.Interpolations: AbstractInterpolation
import ImmersedLayers.ConstrainedSystems.RecursiveArrayTools: ArrayPartition 

const DEFAULT_DT = 0.01
const DEFAULT_DT_STREAK = 0.01
const DEFAULT_ALG = Euler()
const DEFAULT_T_DURATION = 3.0

"""
    Trajectories

Type returned by trajectory calculations, containing a set of one or more trajectories.
For an instance `traj` of this type, the
trajectory time array can be returned with `traj.t`. The number of trajectories is returned
by `traj.np`. Any trajectory contained in the set can be obtained with `traj[p]`, where
`p` must be `0 < p <= traj.np`. 
"""
struct Trajectories{TT,XT}
  np :: Integer
  t :: TT
  xhistory :: Vector{XT}
  yhistory :: Vector{XT}
end

function Trajectories(sol::ODESolution)
  xhist, yhist = _trajectories(sol.u)
  Trajectories(length(xhist[1]),sol.t,xhist,yhist)
end

function _trajectories(u::Vector{T}) where {T<:ArrayPartition}
  xhist = map(s -> s.x[1],u)
  yhist = map(s -> s.x[2],u)
  return xhist, yhist
end

function _trajectories(u::Vector{T}) where {T<:Vector}
  xhist = map(s -> s[1],u)
  yhist = map(s -> s[2],u)
  return xhist, yhist
end

Base.size(traj::Trajectories) = traj.np

Base.length(traj::Trajectories) = traj.np

Base.getindex(traj::Trajectories,k::Integer) = _pick_trajectory(traj,k)

function _pick_trajectory(traj::Trajectories,pnum::Integer)
  @assert pnum > 0 && pnum <= traj.np "Unavailable trajectory number"
  return map(x -> x[pnum],traj.xhistory), map(y -> y[pnum],traj.yhistory)
end

## APIs ##



"""
    displacement_field(vr::Vector{Tuple{AbstractInterpolation,AbstractInterpolation}},tr::AbstractVector,x0::ScalarGridData,y0::ScalarGridData,Trange::Tuple[;Δt=step(tr),alg=Euler()])

Calculate the displacement of particles initally at coordinates `x0` and `y0` over the range of times `Trange = (ti,tf)`, using the vector of spatially-interpolated
velocity fields in `vr` (with corresponding times `tr`). The final time in `Trange` can be earlier than the initial time if backward trajectories are desired.
The optional keyword arguments are `Δt`, the time step size (which defaults to the step size in `tr`, but could be an integer multiple larger than 1). 
"""
function displacement_field(vr::Vector{Tuple{T,T}},tr::StepRangeLen,x0::ScalarGridData,y0::ScalarGridData,Trange::Tuple;Δt::Real=step(tr),alg=ILMPostProcessing.DEFAULT_ALG,kwargs...) where T<:AbstractInterpolation
  ti, tf = Trange
  traj = compute_trajectory(vr,tr,(x0,y0),Trange,alg=Euler(),saveat=[tf])

  xf, yf = traj.xhistory[1], traj.yhistory[1]

  return xf, yf
end


"""
   compute_trajectory(vr::Vector{Tuple{AbstractInterpolation,AbstractInterpolation}},tr::AbstractVector,X₀,Trange::Tuple[;Δt=step(tr),alg=Euler()])

Calculate the trajectories of particles with initial location(s) `X₀`. The argument
`vr` is a vector of spatially-interpolated velocity fields, `tr` is the corresponding time array,
`X₀` can be specified as either a single vector `[x0,y0]`, a vector of vectors specifying
x, y pairs, or a tuple of vectors or arrays specifying x and y positions, respectively,
for multiple tracer particles. `Trange` is a tuple of the starting
and ending integration times. The optional keyword arguments are `Δt`, the time step size (which defaults to the step size in `tr`, but could be an integer multiple larger than 1). The output is the solution
structure for the `OrdinaryDiffEq` package.
"""
function compute_trajectory(vr::Vector{Tuple{T,T}},tr::StepRangeLen,X₀,Trange::Tuple;Δt::Real=step(tr),alg=DEFAULT_ALG,kwargs...) where T<:AbstractInterpolation
  @assert length(vr) == length(tr) "Supplied time array must be same length as supplied velocity array"
  ti, tf = _check_times(tr,Trange,Δt)
  _dt = sign(tf-ti)*abs(Δt)

  u0 = _prepare_initial_conditions(X₀)
  vfcn!(dR,R,p,t) = _vfcn_interpolated_series!(dR,R,p,t,vr,tr)

  sol = _solve_trajectory(vfcn!,u0,Trange,_dt,alg;kwargs...)

  return Trajectories(sol)
  
end



"""
   compute_trajectory(u,v,X₀,Trange::Tuple[;Δt=0.01])

Calculate the trajectory of a tracer particle with initial location(s) `X₀`, which
can be specified as either a single vector `[x0,y0]`, a vector of vectors specifying
x, y pairs, or a tuple of vectors or arrays specifying x and y positions, respectively,
for multiple tracer particles. The arguments
`u` and `v` are either interpolated velocity field components from a computational solution
or are functions. If they are functions, then each of them should be of the form `u(x,y,t)`
and `v(x,y,t)`; `Trange` is a tuple of the initial and final time of integration (and the final time
can be earlier than the initial time if backward trajectories are desired); and `Δt` is the
time step size, which defaults to 0.001. The output is the solution
structure for the `OrdinaryDiffEq` package (or, for multiple particles, a vector
of such solution structures).
"""
function compute_trajectory(ufield::T,
                            vfield::T,
                            X₀,Trange::Tuple;Δt::Real=DEFAULT_DT,alg=DEFAULT_ALG,kwargs...) where {T<:AbstractInterpolation}

  u0 = _prepare_initial_conditions(X₀)
  sol = _compute_trajectory(ufield,vfield,u0,Trange,Δt,alg;kwargs...)
  return Trajectories(sol)

end

function compute_trajectory(ufcn::Function,
                            vfcn::Function,
                            X₀,Trange::Tuple;Δt::Real=DEFAULT_DT,alg=DEFAULT_ALG,kwargs...)

  ti, tf = Trange
  _dt = sign(tf-ti)*abs(Δt)
  velfcn(R,p,t) = _is_autonomous_velocity(ufcn) ? _vfcn_autonomous(R,p,t,ufcn,vfcn) : _vfcn_nonautonomous(R,p,t,ufcn,vfcn)

  u0 = _prepare_initial_conditions(X₀)
  sol = _solve_trajectory(velfcn,u0,Trange,_dt,alg;kwargs...)

  return Trajectories(sol)

end



#######

function _compute_trajectory(vr::Vector{Tuple{T,T}},tr::StepRangeLen,X₀,Trange::Tuple,Δt,alg;kwargs...) where T<:AbstractInterpolation
  @assert length(vr) == length(tr) "Supplied time array must be same length as supplied velocity array"
  ti, tf = _check_times(tr,Trange,Δt)
  _dt = sign(tf-ti)*abs(Δt)

  vfcn!(dR,R,p,t) = _vfcn_interpolated_series!(dR,R,p,t,vr,tr)

  sol = _solve_trajectory(vfcn!,X₀,Trange,_dt,alg;kwargs...)
  return sol
  
end

function _compute_trajectory(ufield::AbstractInterpolation{T,2},vfield::AbstractInterpolation{T,2},
                              u0,Trange,Δt,alg;kwargs...) where {T}

  ti, tf = Trange
  _dt = sign(tf-ti)*abs(Δt)

  vfcn!(dR,R,p,t) = _vfcn_autonomous!(dR,R,p,t,ufield,vfield)

  sol = _solve_trajectory(vfcn!,u0,Trange,_dt,alg;kwargs...)
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
                            X₀::Vector{S},srange::Tuple,t::Real;Δs::Real=DEFAULT_DT,alg=DEFAULT_ALG,kwargs...) where {S<:Real}

  velfcn(R,p,s) = _vfcn_nonautonomous_frozentime(R,p,s,ufcn,vfcn)

  sol = _solve_streamline(velfcn,X₀,srange,Δs,t,alg;kwargs...)
  return sol

end

function compute_streamline(ufield,vfield,
   pts::Vector{Vector{S}},srange::Tuple,t::Real;Δs=DEFAULT_DT,alg=DEFAULT_ALG,kwargs...) where {S<:Real}

  sol_array = ODESolution[]
  for X₀ in pts
    sol = compute_streamline(ufield,vfield,X₀,srange,t;Δs=Δs,alg=alg,kwargs...)
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
function compute_streakline(u,v,X₀::Vector{S},t;τmin = t-DEFAULT_T_DURATION, Δtstreak::Real=DEFAULT_DT_STREAK, Δttraj::Real=DEFAULT_DT,alg=DEFAULT_ALG,kwargs...) where {S<:Real}
  τstreak = τmin:Δtstreak:t
  xstreak = zeros(length(τstreak))
  ystreak = zeros(length(τstreak))

  for (i,τ) in enumerate(τstreak)
    traj = compute_trajectory(u,v,X₀,(τ,t);Δt = Δttraj,alg=alg,kwargs...)
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

function _solve_trajectory(vfcn,u0,Trange,Δt,alg; kwargs...)
  Path = ODEProblem(vfcn,u0,Trange)
  sol = OrdinaryDiffEq.solve(Path,alg; dt = Δt, maxiters = 1e8, adaptive = false, dense = false, kwargs...)
end

function _solve_streamline(vfcn,u0,Trange,Δt,p,alg; kwargs...)
  Path = ODEProblem(vfcn,u0,Trange,p)
  sol = OrdinaryDiffEq.solve(Path,alg; dt = Δt, maxiters = 1e8, adaptive = false, dense = false, kwargs...)
end

## Right-hand side functions for trajectories ##

function _vfcn_autonomous!(dR::ArrayPartition,R::ArrayPartition,p,t,u,v)
 dR.x[1] .= u.(R.x[1],R.x[2])
 dR.x[2] .= v.(R.x[1],R.x[2])

 return dR
end

function _vfcn_autonomous!(dR,R,p,t,u,v)
  dR[1] = u(R[1],R[2])
  dR[2] = v(R[1],R[2])
  return dR
 end

 function _vfcn_autonomous(R::ArrayPartition,p,t,u,v)
  dR = similar(R)
  dR.x[1] .= u.(R.x[1],R.x[2])
  dR.x[2] .= v.(R.x[1],R.x[2])
  return dR
 end

 function _vfcn_autonomous(R,p,t,u,v)
  dR = similar(R)
  dR[1] = u(R[1],R[2])
  dR[2] = v(R[1],R[2])
 
  return dR
 end

function _vfcn_nonautonomous(R::ArrayPartition,p,t,u,v)
 dR = similar(R)
 dR.x[1] .= u.(R.x[1],R.x[2],Ref(t))
 dR.x[2] .= v.(R.x[1],R.x[2],Ref(t))

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

function _vfcn_interpolated_series!(dR,R,p,t,vr,tr)
  jr = searchsortedfirst(tr,t)
  u, v = vr[jr]
  dR.x[1] .= u.(R.x[1],R.x[2])
  dR.x[2] .= v.(R.x[1],R.x[2])
  return dR
end

## For computing fields along trajectories ##

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
   

function _check_times(tr,Trange,Δt)
    ti, tf = Trange
    ji, jf = searchsortedfirst(tr,ti), searchsortedfirst(tr,tf)
    @assert tr[ji] ≈ ti "First entry in time range is not in supplied time array"
    @assert tr[jf] ≈ tf "Last entry in time range is not in supplied time array"
    @assert abs(Δt/step(tr)) >= 1 "Supplied time step size must be >= than step size in supplied time array"
    @assert mod(Δt,step(tr)) ≈ 0 "Supplied time step size must be integer multiple of step size in supplied time array"
    return ti, tf
end
       

function _prepare_initial_conditions(X₀::Tuple)
  x0, y0 = X₀
  return ArrayPartition(deepcopy(x0),deepcopy(y0))
end

function _prepare_initial_conditions(X₀::Vector{Vector{S}}) where {S<:Real}
  x0 = map(pt -> pt[1],X₀)
  y0 = map(pt -> pt[2],X₀)
  return ArrayPartition(x0,y0)
end

function _prepare_initial_conditions(X₀::Vector{S}) where {S<:Real}
  return X₀
end

function _is_autonomous_velocity(u::Function)
  m = first(methods(u))
  # autonomous will have three (fcn name, x, y), non-autonomous will have four (time, as well)
  is_aut = m.nargs == 3 ? true : false
  return is_aut
end