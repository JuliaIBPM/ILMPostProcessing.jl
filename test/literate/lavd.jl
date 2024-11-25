#=
# Lagrangian-averaged vorticity deviation (LAVD)
In this example, we will compute the Lagrangian-averaged vorticity deviation (LAVD) field for a flat plate undergoing a 45-degree pitch-up maneuver. The original simluation of the flat plate was used by Wang and Eldredge 2012 (https://doi.org/10.1007/s00162-012-0279-5). The results are compared with Huang and Green 2016 (https://arc.aiaa.org/doi/10.2514/6.2016-2082). The theory behind LAVD can be found in Haller et al. 2016 (https://doi.org/10.1017/jfm.2016.151). A MATLAB package for computing LAVD and extracting coherent vortices can be found here (https://github.com/Hadjighasem/Lagrangian-Averaged-Vorticity-Deviation-LAVD).
=#

#md # ```@meta
#md # CurrentModule = ILMPostProcessing
#md # ```

using ILMPostProcessing
using ViscousFlow
using Plots
using Statistics

#=
# Viscous Flow of Pitching Flat Plate
## Problem Specification
For faster compututation and testing purposes, the Reynolds number is set to 100 as opposed to 1000 in Huang's paper. The grid Re is also set to 4.0. If better resolution is desired, try grid Re = 3.0. The domain of interest is from x = -0.5 to x = 5.5, but it is set from x = -3.0 to x = 5.5 since the velocity and vorticity fields ahead of the flat plate are required to compute LAVD.
=#

my_params = Dict() 
my_params["Re"] = 100
my_params["grid Re"] = 4.0
my_params["freestream speed"] = 1.0
my_params["freestream angle"] = 0.0

xlim = (-3.0,5.5)
ylim = (-2.0,1.0)
g = setup_grid(xlim,ylim,my_params)
Δs = surface_point_spacing(g,my_params)

#=
## Set up Body
Create a rectangle of length 1.0 and thickness 0.023.
=#

Lp = 1.0
body = Rectangle(Lp/2,0.023/2,Δs)
bl = BodyList([body])

#=
## Set the Body Motion
Create a smooth position ramp of 45 degreess of the flat plate's angle of attack about its leading edge.
=#

vel = 45pi/180  ## nominal ramp velocity
Δx = -45pi/180 ## change in position
t0 = 1.0 ## time of ramp start
k = SmoothRampDOF(vel,Δx,t0)

#=
Plot the ramp. 
=#

t = range(0,3,length=301)
plot(t,dof_position.(k.(t)),xlims=(0,Inf),label="x")

#=
Create the joint. 
=#

parent_body, child_body = 0, 1
Xp = MotionTransform([0,0],0) # location of joint in inertial system
xpiv = [-0.5,0] # place center of motion at LE
Xc = MotionTransform(xpiv,0)
joint1 = Joint(RevoluteJoint,parent_body,Xp,child_body,Xc,[k])
m = RigidBodyMotion([joint1],bl)

x = init_motion_state(bl,m)
update_body!(bl,x,m)
plot(bl,xlim=xlim,ylim=ylim)

#=
Animate the motion
=#

macro animate_motion(b,m,dt,tmax,xlim,ylim)
    return esc(quote
            bc = deepcopy($b)
            t0, x0 = 0.0, init_motion_state(bc,$m)
            dxdt = zero(x0)
            x = copy(x0)

            @gif for t in t0:$dt:t0+$tmax
                motion_rhs!(dxdt,x,($m,bc),t)
                global x += dxdt*$dt
                update_body!(bc,x,$m)
                plot(bc,xlim=$xlim,ylim=$ylim)
            end every 5
        end)
end

@animate_motion bl m 0.01 4 (-0.5, 5.5) ylim

#=
## Define the Boundary Condition Functions
=#

function my_vsplus(t,x,base_cache,phys_params,motions)
  vsplus = zeros_surface(base_cache)
  surface_velocity!(vsplus,x,base_cache,motions,t)
  return vsplus
end

function my_vsminus(t,x,base_cache,phys_params,motions)
  vsminus = zeros_surface(base_cache)
  surface_velocity!(vsminus,x,base_cache,motions,t)
  return vsminus
end

bcdict = Dict("exterior" => my_vsplus, "interior" => my_vsminus)

#=
## Construct the system structure
=#

sys = viscousflow_system(g,bl,phys_params=my_params,motions=m,bc=bcdict);
u0 = init_sol(sys)

#=
Check the effective Reynolds number.
=#

Umax, imax, tmax, bmax = maxvelocity(u0,sys)
L = Lp
Re_eff = my_params["Re"]*Umax*L

#=
Initialize the solver.
=#

tspan = (0.0,10.0)
integrator = init(u0,tspan,sys,alg=LiskaIFHERK(saddlesolver=CG))

step!(integrator)
@time step!(integrator,9.0)

#=
Plot the solutions.
=#

sol = integrator.sol
plt = plot(layout = (5,2), size = (1000, 1200), legend=:false)
tsnap = 0.0:1.0:9.0
for (i, t) in enumerate(tsnap)
    plot!(plt[i],vorticity(sol,sys,t),sys,layers=false,title="t = $(round(t,digits=2))",clim=(-5,5),levels=range(-5,5,length=30),color = :RdBu, xlim = (-0.5, 5.5))
    plot!(plt[i],surfaces(sol,sys,t))
end
plt

#=
# Compute LAVD 
=#

#=
## Generate a Sequence of Velocity and Vorticity Fields
This step obtains the computed velocity and vorticity fields at a sequence of times, and stores them as a sequence of interpolatable
fields. This will greatly speed up how we compute the flow properties (i.e. vorticity) along trajectories.
=#

t_start = 0.0
t_end = 8.5
dt = timestep(u0,sys)
tr = t_start:dt:t_end

velxy = velocity_xy(sol,sys,tr) # Vector of interpolatable velocities
velseq = VectorFieldSequence(tr,velxy); # Bundle together with the time array
vortxy = vorticity_xy(sol,sys,tr) 
vortseq = ScalarFieldSequence(tr,vortxy);

#=
## Generate Initial Conditions
Here, we generate a grid of initial locations from which to integrate trajectories.
=#

X_MIN = -0.5
X_MAX = 5.5
Y_MIN = -2.0
Y_MAX = 1.0
dx = 0.04
lavdgrid = PhysicalGrid((X_MIN,X_MAX),(Y_MIN,Y_MAX),dx)
lavd_cache = SurfaceScalarCache(lavdgrid)
x0, y0 = x_grid(lavd_cache), y_grid(lavd_cache)

#=
## Solve the IVP and Generate LAVD Fields
=#

#=
Compute trajectories
=#

T = -2.0
t0 = 8.5
t1 = t0 + T

@time traj = compute_trajectory(velseq, (x0, y0), (t0, t1));

#=
Compute LAVD 
=#

LAVD = similar(x0)
compute_LAVD!(LAVD, traj, X_MIN, Y_MIN, X_MAX, Y_MAX, vortseq)
plot(LAVD, lavd_cache, colorbar = true, levels = 20)
plot!(surfaces(sol,sys,t0))
#savefig("lavd")

#=
Compute IVD 
=#

IVD = similar(x0)
compute_IVD!(IVD, traj, X_MIN, Y_MIN, X_MAX, Y_MAX, vortseq)
plot(IVD, lavd_cache, colorbar = true, levels = 20)
plot!(surfaces(sol,sys,t0))
#savefig("ivd")

#=
Plot the vorticity fields
=#

plot(vorticity(sol, sys, t0), sys, layers = false,clim = (-5,5),levels = range(-5,5,length=30), colorbar = true, xlim = (-0.5, 5.5))
plot!(surfaces(sol,sys,t0))
#savefig("vorticity")
