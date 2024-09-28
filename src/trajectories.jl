# TRAJECTORY CALCULATION #

import Interpolations: AbstractInterpolation
import RecursiveArrayTools: ArrayPartition 

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

### field sequences ###

abstract type AbstractFieldSequence{FT,TT} end

"""
    VectorFieldSequence(tseq::AbstractVector,vseq)

This type bundles a time vector with a vector of tuples of interpolatable
fields (i.e., each member of the tuple is of type `AbstractInterpolation`
with two spatial coordinate arguments). It is used in trajectory computations
and for plotting fields along trajectories. 
"""
struct VectorFieldSequence{FT,TT} <: AbstractFieldSequence{FT,TT}
  tseq :: TT
  fseq :: Vector{Tuple{FT,FT}}
end

"""
    ScalarFieldSequence(tseq::AbstractVector,sseq)

This type bundles a time vector with a vector of interpolatable
scalar fields (i.e., each element is of type `AbstractInterpolation`
with two spatial coordinate arguments). It is used for plotting fields along
trajectories. 
"""
struct ScalarFieldSequence{FT,TT} <: AbstractFieldSequence{FT,TT}
  tseq :: TT
  fseq :: Vector{FT}
end

Base.step(v::AbstractFieldSequence{F,T}) where {F,T <: AbstractRange} = step(v.tseq)
Base.step(v::AbstractFieldSequence{F,T}) where {F,T <: Vector} = v.tseq[2]-v.tseq[1]


# These functions look up the field in the sequence at the time closest to t
function _instantaneous_vector_field_in_series(v::VectorFieldSequence,t)
  jr = searchsorted(v.tseq,t)
  j1, j2 = first(jr), last(jr)
  jt = abs(v.tseq[j1] - t) <= abs(v.tseq[j2] - t) ? j1 : j2
  return v.fseq[jt]
end

function _instantaneous_scalar_field_in_series(s::ScalarFieldSequence,t)
  jr = searchsorted(s.tseq,t)
  j1, j2 = first(jr), last(jr)
  jt = abs(s.tseq[j1] - t) <= abs(s.tseq[j2] - t) ? j1 : j2
  return s.fseq[jt]
end


## APIs ##



"""
    displacement_field(v::VectorFieldSequence,x0,y0,Trange::Tuple[;Δt=step(tr),alg=Euler()])

Calculate the displacement of particles initally at coordinates `x0` and `y0` over the range of times `Trange = (ti,tf)`, using the sequence of spatially-interpolated
velocity fields in `v`. (One can also provide the vector of velocity fields and the time array as separate arguments).
The final time in `Trange` can be earlier than the initial time if backward trajectories are desired.
The optional keyword arguments are `Δt`, the time step size (which defaults to the step size in `tr`, but could be an integer multiple larger than 1). 
"""
function displacement_field(v::VectorFieldSequence,x0,y0,Trange::Tuple;Δt::Real=step(v),alg=ILMPostProcessing.DEFAULT_ALG,kwargs...)
  ti, tf = Trange
  traj = compute_trajectory(v,(x0,y0),Trange,alg=Euler(),saveat=[tf])

  xf, yf = traj.xhistory[1], traj.yhistory[1]

  return xf, yf
end

displacement_field(vr::Vector{Tuple{T,T}},tr::StepRangeLen,x0,y0,Trange;kwargs...) where T<:AbstractInterpolation = 
        displacement_field(VectorFieldSequence(tr,vr),x0,y0,Trange;kwargs...)


"""
   compute_trajectory(v::VectorFieldSequence,X₀,Trange::Tuple[;Δt=step(tr),alg=Euler()])

Calculate the trajectories of particles with initial location(s) `X₀`. The argument
`v` contains a sequence of spatially-interpolated velocity fields and an associated time array.
(One can also provide the vector of velocity fields and the time array as separate arguments).

`X₀` can be specified as either a single vector `[x0,y0]`, a vector of vectors specifying
x, y pairs, or a tuple of vectors or arrays specifying x and y positions, respectively,
for multiple tracer particles. `Trange` is a tuple of the starting
and ending integration times. The optional keyword arguments are `Δt`, the time step size (which defaults to the step size in `tr`, but could be an integer multiple larger than 1). The output is the solution
structure for the `OrdinaryDiffEq` package.
"""
function compute_trajectory(v::VectorFieldSequence,X₀,Trange::Tuple;Δt::Real=step(v),alg=DEFAULT_ALG,kwargs...)
  ti, tf = _check_times(v.tseq,Trange,Δt)
  tsign = sign(tf-ti)
  _dt = tsign != 0 ? tsign*abs(Δt) : Δt

  u0 = _prepare_initial_conditions(X₀)
  vfcn!(dR,R,p,t) = _vfcn_interpolated_series!(dR,R,p,t,v)

  sol = _solve_trajectory(vfcn!,u0,Trange,_dt,alg;kwargs...)

  return Trajectories(sol)
  
end

compute_trajectory(vr::Vector{Tuple{T,T}},tr::StepRangeLen,X₀,Trange::Tuple;kwargs...) where {T<:AbstractInterpolation} = compute_trajectory(VectorFieldSequence(tr,vr),X₀,Trange;kwargs...)



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
  tsign = sign(tf-ti)
  _dt = tsign != 0 ? tsign*abs(Δt) : Δt

  velfcn(R,p,t) = _is_autonomous_velocity(ufcn) ? _vfcn_autonomous(R,p,t,ufcn,vfcn) : _vfcn_nonautonomous(R,p,t,ufcn,vfcn)

  u0 = _prepare_initial_conditions(X₀)
  sol = _solve_trajectory(velfcn,u0,Trange,_dt,alg;kwargs...)

  return Trajectories(sol)

end



#######

#=
function _compute_trajectory(vr::Vector{Tuple{T,T}},tr::StepRangeLen,X₀,Trange::Tuple,Δt,alg;kwargs...) where T<:AbstractInterpolation
  @assert length(vr) == length(tr) "Supplied time array must be same length as supplied velocity array"
  tsign = sign(tf-ti)
  _dt = tsign != 0 ? tsign*abs(Δt) : Δt

  vfcn!(dR,R,p,t) = _vfcn_interpolated_series!(dR,R,p,t,vr,tr)

  sol = _solve_trajectory(vfcn!,X₀,Trange,_dt,alg;kwargs...)
  return sol
  
end
=#

function _compute_trajectory(ufield::AbstractInterpolation{T,2},vfield::AbstractInterpolation{T,2},
                              u0,Trange,Δt,alg;kwargs...) where {T}

  ti, tf = Trange
  tsign = sign(tf-ti)
  _dt = tsign != 0 ? tsign*abs(Δt) : Δt

  vfcn!(dR,R,p,t) = _vfcn_autonomous!(dR,R,p,t,ufield,vfield)

  sol = _solve_trajectory(vfcn!,u0,Trange,_dt,alg;kwargs...)
  return sol

end

"""
   compute_streamline(u,v,X₀,srange::Tuple,t::Real,[;Δt=0.01])

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
    xtraj, ytraj = traj[1]
    xstreak[i], ystreak[i] = xtraj[end],ytraj[end]
  end
  return Trajectories(1,τstreak,xstreak,ystreak)
end


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

function _field_along_trajectory(sseq::VectorFieldSequence,traj::Trajectories,p,::Val{0})
  xh, yh = traj[p]
  varray = map((x,y,t) -> (vel = _instantaneous_vector_field_in_series(sseq,t); tuple(vel[1](x,y), vel[2](x,y))),xh,yh,traj.t)
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