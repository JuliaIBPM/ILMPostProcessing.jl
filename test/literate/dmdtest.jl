#=
# Dynamic mode decomposition (DMD)
In this example, we will demonstrate the use of dynamic mode decomposition (DMD) for
decomposing a simple linear dynamical system.

This example is inspired from example 1 of

M.S. Hemati, C.W. Rowley, E.A. Deem, and L.N. Cattafesta
   ``De-biasing the dynamic mode decomposition for 
     applied Koopman spectral analysis of noisy datasets,''
     Theoretical and Computational Fluid Dynamics (2017).

The example considers a low-rank linear system with 
two undamped modes and one dampled mode. The snapshots taken from the
solution of the linear system are noised up with zero-mean Gaussian noise.
=#

#md # ```@meta
#md # CurrentModule = ILMPostProcessing
#md # ```

using ILMPostProcessing

using LinearAlgebra
using Random
using OrdinaryDiffEq
#!jl using Plots

m = 100 # number of snapshots
n = 250 # number of states
r = 6 # rank of DMD
dt = 0.01 # snapshot time step
meas_cov = 0.05 # measurement noise covariance
init_cov = 0.1 # initial condition covariance

#=
Specify characteristic frequencies and growth/decay rates
associated with continuous-time dynamics.
The DMD rank should be set equal to twice the number of modes
(since each mode consists of conjugate pairs)
=#
f = [1.0, 2.5, 5.5]
g = [0, 0, -0.3]

#=
Create the right hand side matrix for the continuous linear system
=#
k = 2*length(f)
A = zeros(k,k)
for ii in 1:length(f)
    i1, i2 = 2*ii-1, 2*ii
    Ai = view(A,i1:i2,i1:i2)
    Ai .= [g[ii] 2π*f[ii]; -2π*f[ii] g[ii]]
end

#=
The true eigenvalues of the system 
=#
true_evals = exp.(eigvals(A)*dt)

#=
Right-hand side of linear system of equations
=#
dynsys(x,p,t) = A*x

#=
### Solve the linear system
Set up a random initial condition with elements drawn from N(1,init_cov)
and solve the problem.
=#
x0 = 1 .+ randn(k)*sqrt(init_cov)

tspan = (0,dt*m)
prob = ODEProblem(dynsys,x0,tspan)
sol = solve(prob,Tsit5(),saveat=dt);

#=
For DMD, use the solution snapshots, but
randomly rotate them and apply noise to each
=#
Q, _ = qr(randn(n,k))
getsnaps(x) = Q*x .+ sqrt(meas_cov)*randn(n)
snaps = map(x -> getsnaps(x),sol.u);


#=
Now perform DMD
=#
dmdmodes = dmd(snaps,r)

#!jl scatter(real(true_evals),imag(true_evals),ratio=1,xlim = (0.7,1.1),ylim=(0,0.4))
#!jl scatter!(real(dmdmodes.evals),imag(dmdmodes.evals))
#!jl θ = range(0,2π,length=100);
#!jl plot!(cos.(θ),sin.(θ))

#=
### Compare the true and DMD-computed eigenvalues
Note that these may not be ordered the same, so we have to
also determine how to permute the order of them to compare
corresponding eigenvalues. We then compute the l2 error
=#
vals, idex = findmin(abs2.(true_evals .- transpose(dmdmodes.evals)),dims=2)
err = sqrt(sum(vals))

#jl @test err < 0.05

#md # ## DMD functions
#md # ```@docs
#md # dmd
#md # ```
