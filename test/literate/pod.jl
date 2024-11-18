#=
# Proper orthogonal decomposition (POD)
In this example, we will demonstrate the use of proper orthogonal decomposition (POD) for
decomposing a flow field into basis modes.
=#

#md # ```@meta
#md # CurrentModule = ILMPostProcessing
#md # ```

using ILMPostProcessing
#!jl using Plots


#=
## Get the flow field data
First, we need some flow field data to analyze. For this purpose, we will use [ViscousFlow.jl](https://github.com/JuliaIBPM/ViscousFlow.jl)
to get snapshots of the flow for a flat plate at 30 degrees angle of attack at
Reynolds number 100.
=#
using ViscousFlow

my_params = Dict()
my_params["Re"] = 100
my_params["freestream speed"] = 1.0 # in x-dir
my_params["freestream angle"] = 0.0 # relative to horizontal

xlim = (-1.0, 5.0)
ylim = (-2.0, 2.0)
my_params["grid Re"] = 4.0

g = setup_grid(xlim, ylim, my_params)

Δs = surface_point_spacing(g, my_params)
body = Plate(1.0, Δs)

cent = [0.5,0.0] 
α = -30π/180 
X = MotionTransform(cent, α)
joint = Joint(X)

m = RigidBodyMotion(joint, body)
x = init_motion_state(body, m) 
update_body!(body, x, m)
sys = viscousflow_system(g, body, phys_params = my_params, motions = m);

u0 = init_sol(sys)
tspan = (0.0, 20.0)
integrator = init(u0, tspan, sys)

## Solve to 10 convective time units
step!(integrator, 10)

#=
## Assemble snapshots of the velocity field from the solution data
Here, we make use of the capability of the `velocity` function to
generate an array of velocity fields at a range of times. We will
save every 5th time step in this array.
=#
sol = integrator.sol
tpod = sol.t[2:5:end]
X = velocity(sol, sys, tpod);

#=
## Perform the POD
The POD is simply performed with the `PODModes` function. This provides
a structure containing the modes (`phi`), the expansion coefficients (`a`), and the modal
energies (`lambda`). By default, `PODModes` retains 99% of the energy. This can be changed
with the optional argument `tolerance`.
=#
modes = pod(X);

#=
 The `a` array is of size $N_t \times r$, where $N_t$ is the number of time values,
and $r$ is the number of modes. The modes are ordered from highest energy to lowest energy.
=#
modes.a

#=
In this case, 7 modes were retained, at 51 times.
=#

#=
If we wanted to re-assemble the modes and coefficients to recover the flow at some time instant, we could
use the `mapreduce` function, e.g.,
=#
vel_assemble = mapreduce((aj, phi_j) -> aj .* phi_j, +, modes.a[end,:], modes.phi) + modes.Xmean

#=
In this last line, `modes.a[end,:]` obtains the expansion coefficients at the last time
available.
=#

#=
Let's print the first mode, and the corresponding history of the modal coefficient in the decomposition
=#
#!jl plot(layout=[2;1],plot(modes.phi[1].u,sys,title="u"),
#!jl    plot(modes.phi[1].v,sys,title="v"),
#!jl    plot(tpod,modes.a[:,1],xlim=(0,Inf),xlabel="\$t\$",ylabel="\$a_1(t)\$"))

#=
The energy associated with this mode is
=#
modes.lambda[1]

#=
Now let's print the $r$th mode, and the history of the coefficient in the decomposition
=#
#!jl plot(layout=[2;1],plot(modes.phi[end].u,sys,title="u"),
#!jl    plot(modes.phi[end].v,sys,title="v"),
#!jl    plot(tpod,modes.a[:,end],xlim=(0,Inf),xlabel="\$t\$",ylabel="\$a_r(t)\$"))

#=
The energy associated with this mode is
=#
modes.lambda[end]