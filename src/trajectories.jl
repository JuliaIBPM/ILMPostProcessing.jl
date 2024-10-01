# TRAJECTORY CALCULATION #

import Interpolations: AbstractInterpolation
import RecursiveArrayTools: ArrayPartition 

const DEFAULT_DT = 0.01
const DEFAULT_DT_STREAK = 0.01
const DEFAULT_ALG = Euler()
const DEFAULT_ADAPTIVE_ALG = Tsit5()
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

### field sequences ###

abstract type AbstractFieldSequence{FT,TT} end

"""
    VectorFieldSequence(t::AbstractVector,v)

This type bundles a time vector `t` with a vector `v` of tuples of interpolatable
fields (i.e., each member of the tuple is of type `AbstractInterpolation`
with two spatial coordinate arguments). It is used in trajectory computations
and for plotting fields along trajectories. 
"""
struct VectorFieldSequence{FT,TT} <: AbstractFieldSequence{FT,TT}
  t :: TT
  v :: Vector{Tuple{FT,FT}}
end

"""
    ScalarFieldSequence(t::AbstractVector,s)

This type bundles a time vector `t` with a vector `s` of interpolatable
scalar fields (i.e., each element is of type `AbstractInterpolation`
with two spatial coordinate arguments). It is used for plotting fields along
trajectories. 
"""
struct ScalarFieldSequence{FT,TT} <: AbstractFieldSequence{FT,TT}
  t :: TT
  v :: Vector{FT}
end

Base.step(seq::AbstractFieldSequence{F,T}) where {F,T <: AbstractRange} = step(seq.t)
Base.step(seq::AbstractFieldSequence{F,T}) where {F,T <: Vector} = seq.t[2]-seq.t[1]


# These functions look up the field in the sequence at the time closest to t
function _instantaneous_vector_field_in_series(seq::VectorFieldSequence,t)
  jr = searchsorted(seq.t,t)
  j1, j2 = first(jr), last(jr)
  jt = abs(seq.t[j1] - t) <= abs(seq.t[j2] - t) ? j1 : j2
  return seq.v[jt]
end

function _instantaneous_scalar_field_in_series(seq::ScalarFieldSequence,t)
  jr = searchsorted(seq.t,t)
  j1, j2 = first(jr), last(jr)
  jt = abs(seq.t[j1] - t) <= abs(seq.t[j2] - t) ? j1 : j2
  return seq.v[jt]
end


## APIs ##



"""
    displacement_field(v::VectorFieldSequence,x0,y0,Trange::Tuple[;dt=step(tr),alg=Euler()])

Calculate the displacement of particles initally at coordinates `x0` and `y0` over the range of times `Trange = (ti,tf)`, using the sequence of spatially-interpolated
velocity fields in `v`. (One can also provide the vector of velocity fields and the time array as separate arguments).
The final time in `Trange` can be earlier than the initial time if backward trajectories are desired.
The optional keyword arguments are `dt`, the time step size (which defaults to the step size in `tr`, but could be an integer multiple larger than 1). 
"""
function displacement_field(vseq::VectorFieldSequence,x0,y0,Trange::Tuple;dt::Real=step(vseq),alg=DEFAULT_ALG,kwargs...)
  ti, tf = Trange
  traj = compute_trajectory(vseq,(x0,y0),Trange;dt=dt,alg=alg,saveat=[tf],kwargs...)

  xf, yf = traj.xhistory[1], traj.yhistory[1]

  return xf, yf
end



displacement_field(vr::Vector{Tuple{T,T}},tr::StepRangeLen,x0,y0,Trange;kwargs...) where T<:AbstractInterpolation = 
        displacement_field(VectorFieldSequence(tr,vr),x0,y0,Trange;kwargs...)

"""
    displacement_field(u::Function,v::Function,x0,y0,Trange::Tuple[;dt=:auto,alg=Euler()])

Calculate the displacement of particles initally at coordinates `x0` and `y0` over the range of times `Trange = (ti,tf)`, using the 
velocity functions `u` and `v`.  These function can either be autonomous (taking only x and y arguments) or non-autonomous, taking
an additional time argument.
The final time in `Trange` can be earlier than the initial time if backward trajectories are desired.
The optional keyword arguments are `dt`, the time step size (which defaults to `:auto`, for adaptive time marching, but could be specified
to override this). The default time marching algorithm is `Tsit5()`. 
"""
function displacement_field(ufcn::Function,vfcn::Function,x0,y0,Trange::Tuple;dt=:auto,alg=DEFAULT_ADAPTIVE_ALG,kwargs...)
  ti, tf = Trange
  traj = compute_trajectory(ufcn,vfcn,(x0,y0),Trange;dt=dt,alg=alg,saveat=[tf],kwargs...)

  xf, yf = traj.xhistory[1], traj.yhistory[1]

  return xf, yf
end

"""
   compute_trajectory(v::VectorFieldSequence,X₀,Trange::Tuple[;dt=step(tr),alg=Euler()])

Calculate the trajectories of particles with initial location(s) `X₀`. The argument
`v` contains a sequence of spatially-interpolated velocity fields and an associated time array.
(One can also provide the vector of velocity fields and the time array as separate arguments).

`X₀` can be specified as either a single vector `[x0,y0]`, a vector of vectors specifying
x, y pairs, or a tuple of vectors or arrays specifying x and y positions, respectively,
for multiple tracer particles. `Trange` is a tuple of the starting
and ending integration times. The optional keyword arguments are `dt`, the time step size (which defaults to the step size in `tr`, but could be an integer multiple larger than 1). The output is the solution
structure for the `OrdinaryDiffEq` package.
"""
function compute_trajectory(vseq::VectorFieldSequence,X₀,Trange::Tuple;dt::Real=step(vseq),alg=DEFAULT_ALG,kwargs...)
  ti, tf = _check_times(vseq.t,Trange,dt)
  _dt, _autodt = _standardize_time_step(ti,tf,dt)

  u0 = _prepare_initial_conditions(X₀)
  vfcn!(dR,R,p,t) = _vfcn_interpolated_series!(dR,R,p,t,vseq)

  sol = _solve_trajectory(vfcn!,u0,Trange,_dt,alg,Val(false);kwargs...)

  return Trajectories(sol)
  
end

compute_trajectory(vr::Vector{Tuple{T,T}},tr::StepRangeLen,X₀,Trange::Tuple;kwargs...) where {T<:AbstractInterpolation} = compute_trajectory(VectorFieldSequence(tr,vr),X₀,Trange;kwargs...)



"""
   compute_trajectory(u,v,X₀,Trange::Tuple[;dt=0.01])

Calculate the trajectory of a tracer particle with initial location(s) `X₀`, which
can be specified as either a single vector `[x0,y0]`, a vector of vectors specifying
x, y pairs, or a tuple of vectors or arrays specifying x and y positions, respectively,
for multiple tracer particles. The arguments
`u` and `v` are either interpolated velocity field components from a computational solution
or are functions. If they are functions, then each of them should be of the form `u(x,y,t)`
and `v(x,y,t)`; `Trange` is a tuple of the initial and final time of integration (and the final time
can be earlier than the initial time if backward trajectories are desired); and `dt` is the
time step size, which defaults to 0.001. The output is the solution
structure for the `OrdinaryDiffEq` package (or, for multiple particles, a vector
of such solution structures).
"""
function compute_trajectory(ufield::T,
                            vfield::T,
                            X₀,Trange::Tuple;dt=:auto,alg=DEFAULT_ALG,kwargs...) where {T<:AbstractInterpolation}

  _dt, _autodt = _standardize_time_step(Trange...,dt)
  vfcn!(dR,R,p,t) = _vfcn_autonomous!(dR,R,p,t,ufield,vfield)
  u0 = _prepare_initial_conditions(X₀)
  sol = _solve_trajectory(vfcn!,u0,Trange,_dt,alg,Val(_autodt);kwargs...)

  return Trajectories(sol)

end

function compute_trajectory(ufcn::Function,
                            vfcn::Function,
                            X₀,Trange::Tuple;dt=:auto,alg=DEFAULT_ALG,kwargs...)

  
  _dt, _autodt = _standardize_time_step(Trange...,dt)
  #velfcn(R,p,t) = _is_autonomous_velocity(ufcn) ? _vfcn_autonomous(R,p,t,ufcn,vfcn) : _vfcn_nonautonomous(R,p,t,ufcn,vfcn)
  velfcn!(dR,R,p,t) = _is_autonomous_velocity(ufcn) ? _vfcn_autonomous!(dR,R,p,t,ufcn,vfcn) : _vfcn_nonautonomous!(dR,R,p,t,ufcn,vfcn)

  u0 = _prepare_initial_conditions(X₀)
  sol = _solve_trajectory(velfcn!,u0,Trange,_dt,alg,Val(_autodt);kwargs...)

  return Trajectories(sol)

end

compute_trajectory(velfield::Tuple{T,T},a...;kwargs...) where {T<:Union{AbstractInterpolation,Function}} = compute_trajectory(velfield...,a...;kwargs...)


function _standardize_time_step(ti,tf,dt::Real)
  tsign = sign(tf-ti)
  _dt = tsign != 0 ? tsign*abs(dt) : dt
  return _dt, false
end

function _standardize_time_step(ti,tf,dt::Symbol)
  return dt, true
end


#######


"""
   compute_streamline(u,v,X₀,srange::Tuple,t::Real,[;dt=0.01])

Calculate the streamline(s) passing through location(s) `X₀`, which
can be specified as either a single vector `[x0,y0]`, a vector of vectors specifying
x, y pairs, or a tuple of vectors or arrays specifying x and y positions, respectively,
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
                            X₀,srange::Tuple,t::Real;Δs::Real=DEFAULT_DT,alg=DEFAULT_ALG,kwargs...)

  velfcn(R,p,s) = _is_autonomous_velocity(ufcn) ? _vfcn_autonomous(R,p,s,ufcn,vfcn) : _vfcn_nonautonomous_frozentime(R,p,s,ufcn,vfcn)

  u0 = _prepare_initial_conditions(X₀)
  sol = _solve_streamline(velfcn,u0,srange,Δs,t,alg;kwargs...)
  return Trajectories(sol)

end

function compute_streamline(ufield::AbstractInterpolation{T,2},vfield::AbstractInterpolation{T,2},
      X₀,srange::Tuple,t::Real;Δs=DEFAULT_DT,alg=DEFAULT_ALG,kwargs...) where {T}

   vfcn!(dR,R,p,s) = _vfcn_autonomous!(dR,R,p,s,ufield,vfield)
   u0 = _prepare_initial_conditions(X₀)
   sol = _solve_streamline(vfcn!,u0,srange,Δs,t,alg;kwargs...)
   return Trajectories(sol)

end

"""
   compute_streakline(u,v,X₀::Vector,t[;τmin = t-3.0, dtstreak=0.01,dttraj=0.001]) -> Vector, Vector

Calculate a streakline at time `t` for a velocity field `u` and `v`, based on an injection
point `X₀`. The end of the streakline is set by `τmin`, the earliest time
at which a particle passed through the injection point. It defaults to 3 time
units before the current instant `t`. The time step size `dtstreak` sets the resolution
of the streakline (i.e., how often the particles are sampled along the streakline).
It returns arrays of the x and y coordinates of the streakline.
"""
function compute_streakline(vel,X₀::Vector{S},t;τmin = t-DEFAULT_T_DURATION, dtstreak::Real=DEFAULT_DT_STREAK, dttraj::Real=DEFAULT_DT,alg=DEFAULT_ALG,kwargs...) where {S<:Real}
  τstreak = τmin:dtstreak:t
  xstreak = zeros(length(τstreak))
  ystreak = zeros(length(τstreak))

  for (i,τ) in enumerate(τstreak)
    traj = compute_trajectory(vel,X₀,(τ,t);dt = dttraj,alg=alg,kwargs...)
    xtraj, ytraj = traj[1]
    xstreak[i], ystreak[i] = xtraj[end],ytraj[end]
  end
  return Trajectories(1,τstreak,xstreak,ystreak)
end

compute_streakline(u::T,v::T,a...;kwargs...) where {T<:Union{AbstractInterpolation,Function}} =
    compute_streakline((u,v),a...;kwargs...)



"""
  field_along_trajectory(f,traj::Trajectories,p::Integer[,deriv=0])

Evaluate field `f` (given as grid data) along the trajectory number `p` in the 
trajectories specified by `traj`.
The output is the history of `f` along this trajectory. If `f` is a vector field,
then the component histories are output as a tuple. If `deriv=1`, then it
computes the time derivative of the field along the trajectory. The default
is `deriv=0` (no derivative).
"""
field_along_trajectory(d,traj,p;deriv=0) = _field_along_trajectory(d,traj,p,Val(deriv))


## Internal helper functions ##

function _solve_trajectory(vfcn,u0,Trange,dt,alg, ::Val{false}; kwargs...)
  Path = ODEProblem(vfcn,u0,Trange)
  sol = OrdinaryDiffEq.solve(Path,alg; dt = dt, maxiters = 1e8, adaptive = false, dense = false, kwargs...)
end

function _solve_trajectory(vfcn,u0,Trange,dt,alg, ::Val{true}; kwargs...)
  Path = ODEProblem(vfcn,u0,Trange)
  sol = OrdinaryDiffEq.solve(Path,alg; kwargs...)
end

function _solve_streamline(vfcn,u0,Trange,dt,p,alg,::Val{false}; kwargs...)
  Path = ODEProblem(vfcn,u0,Trange,p)
  sol = OrdinaryDiffEq.solve(Path,alg; dt = dt, maxiters = 1e8, adaptive = false, dense = false, kwargs...)
end

function _solve_streamline(vfcn,u0,Trange,dt,p,alg,::Val{true}; kwargs...)
  Path = ODEProblem(vfcn,u0,Trange,p)
  sol = OrdinaryDiffEq.solve(Path,alg; maxiters = 1e8, dense = false, kwargs...)
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

 function _vfcn_nonautonomous!(dR::ArrayPartition,R::ArrayPartition,p,t,u,v)
  dR.x[1] .= u(R.x[1],R.x[2],t)
  dR.x[2] .= v(R.x[1],R.x[2],t)
  return dR
 end

 function _vfcn_nonautonomous!(dR,R,p,t,u,v)
  dR[1] = u(R[1],R[2],t)
  dR[2] = v(R[1],R[2],t)
 
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


 function _vfcn_nonautonomous_frozentime(R::ArrayPartition,p,t,u,v)
  dR = similar(R)
  dR.x[1] .= u.(R.x[1],R.x[2],Ref(p))
  dR.x[2] .= v.(R.x[1],R.x[2],Ref(p))
 
  return dR
 end
 
 function _vfcn_nonautonomous_frozentime(R,p,t,u,v)
   dR = similar(R)
   dR[1] = u(R[1],R[2],p)
   dR[2] = v(R[1],R[2],p)
  
   return dR
  end

function _vfcn_interpolated_series!(dR::ArrayPartition,R::ArrayPartition,p,t,vr)
  u, v = _instantaneous_vector_field_in_series(vr,t)
  dR.x[1] .= u.(R.x[1],R.x[2])
  dR.x[2] .= v.(R.x[1],R.x[2])
  return dR
end

function _vfcn_interpolated_series!(dR,R,p,t,vr)
  u, v = _instantaneous_vector_field_in_series(vr,t)
  dR[1] = u(R[1],R[2])
  dR[2] = v(R[1],R[2])
  return dR
end


## For computing fields along trajectories ##

function _field_along_trajectory(v::Tuple{T,T},traj::Trajectories,p,::Val{0}) where T<:Union{AbstractInterpolation,Function}
    vfield_x, vfield_y = v
   
    xh, yh = traj[p]
    vx_traj = vfield_x.(xh,yh)
    vy_traj = vfield_y.(xh,yh)
    return vx_traj, vy_traj
end
   
function _field_along_trajectory(sfield::T,traj::Trajectories,p,::Val{0}) where T<:Union{AbstractInterpolation,Function}
   
    xh, yh = traj[p]
    s_traj = sfield.(xh,yh)
    return s_traj
end

function _field_along_trajectory(vseq::VectorFieldSequence,traj::Trajectories,p,::Val{0})
  xh, yh = traj[p]
  varray = map((x,y,t) -> (vel = _instantaneous_vector_field_in_series(vseq,t); tuple(vel[1](x,y), vel[2](x,y))),xh,yh,traj.t)
  return  map(v -> v[1],varray), map(v -> v[2],varray)
end

function _field_along_trajectory(sseq::ScalarFieldSequence,traj::Trajectories,p,::Val{0})
    xh, yh = traj[p]
    return  map((x,y,t) -> (f = _instantaneous_scalar_field_in_series(sseq,t); f(x,y)),xh,yh,traj.t)
end

function _field_along_trajectory(v::Tuple{T,T},traj::Trajectories,p,::Val{1}) where T<:Union{AbstractInterpolation,Function}
       utraj, vtraj = _field_along_trajectory(v,traj,p,Val(0))
       return ddt(utraj,traj.t), ddt(vtraj,traj.t)
end
   
_field_along_trajectory(s::T,traj::Trajectories,p,::Val{1}) where T<:Union{AbstractInterpolation,Function,ScalarFieldSequence} =
       ddt(_field_along_trajectory(s,traj,p,Val(0)),traj.t,mydiff=:backward_diff)



_evaluate_field_function(x,y,t,field::Function,::Val{4}) = field(x,y,t)
_evaluate_field_function(x,y,t,field::Function,::Val{3}) = field(x,y)


function _check_times(tr,Trange,dt)
    ti, tf = Trange
    ji, jf = searchsortedfirst(tr,ti), searchsortedfirst(tr,tf)
    @assert tr[ji] ≈ ti "First entry in time range is not in supplied time array"
    @assert tr[jf] ≈ tf "Last entry in time range is not in supplied time array"
    @assert abs(dt/step(tr)) >= 1 "Supplied time step size must be >= than step size in supplied time array"
    @assert mod(dt,step(tr)) ≈ 0 "Supplied time step size must be integer multiple of step size in supplied time array"
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