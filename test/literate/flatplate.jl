using ViscousFlow
# using Plots
using Statistics
using LinearAlgebra
using ILMPostProcessing

# SPECIFY
# define the parameters
my_params = Dict()
my_params["Re"] = 100
my_params["freestream speed"] = 1.0 # in x-dir
my_params["freestream angle"] = 0.0 # relative to horizontal

# DISCRETIZE
# setup grid parameters
xlim = (-1.0, 5.0)
ylim = (-2.0, 2.0)
my_params["grid Re"] = 4.0

# setup grid
g = setup_grid(xlim, ylim, my_params)

# create sys based on grid
sys = viscousflow_system(g, phys_params = my_params)

# define body
Δs = surface_point_spacing(g, my_params)
body = Plate(1.0, Δs) # the plate itself

# transform body
cent = [0.5,0.0] # center of joint wrt inertial sys
α = -30π/180 # angle of joint wrt intertial system (CCW rotation)
X = MotionTransform(cent, α) # in LATEST ViscFlow package (0.6.0)
joint = Joint(X)

m = RigidBodyMotion(joint, body)
x = init_motion_state(body, m) # x is (angle, x coord, y coord) of joint
# set body in place using update_body! command
update_body!(body, x, m)

# plot body to ensure it's popisitioned correctly

sys = viscousflow_system(g, body, phys_params = my_params, motions = m)
# note BODY is now involved in sys and m (RigidBodyMotion)

# initialize solution
u0 = init_sol(sys)

tspan = (0.0, 20.0)
integrator = init(u0, tspan, sys)

# SOLVE, step soln to t = 1.0 convective time units
step!(integrator, 4.5)

sol = integrator.sol

X = createSnapshotData(velocity, sol, sys)
modes = PODModes(X)