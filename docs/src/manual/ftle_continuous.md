```@meta
EditURL = "../../../test/literate/ftle_continuous.jl"
```

# Finite-Time Lyapunov Exponent (FTLE) with Continuous Velocity Field
We present two examples where the FTLE is calculated from known continuous velocity fields as opposed to discrete fields specified at grid points.
The first example demonstrates a simple FTLE analysis and lays the conceptual foundation for later analysis on discrete data. The second example showcases an issue of sliding-window FTLE analysis.

```@meta
CurrentModule = ILMPostProcessing
```

````@example ftle_continuous
using ILMPostProcessing
using ImmersedLayers
using Plots
````

## Example 1: Unsteady Double Gyre
This example replicates the time-dependent double gyre from Shadden 2005. The double gyre velocity field's specification can be found at https://shaddenlab.berkeley.edu/uploads/LCS-tutorial/examples.html#Sec7.1. The results in this simulation is highly similar.
### Generate Initial Conditions for the IVP

````@example ftle_continuous
X_MIN, X_MAX = 0.0, 2.0
Y_MIN, Y_MAX = 0.0, 1.0
dx = 0.002

ftlegrid = PhysicalGrid((X_MIN,X_MAX),(Y_MIN,Y_MAX),dx,optimize=false);
ftle_cache = SurfaceScalarCache(ftlegrid)
x0, y0 = x_grid(ftle_cache), y_grid(ftle_cache)
````

### Define the Double Gyre's Vector Field

````@example ftle_continuous
A = 0.1
epsilon = 0.25
omega = 2 * π / 10

a(s) = epsilon * sin(omega * s)
b(s) = 1 - 2 * epsilon * sin(omega * s)
f(x, t) = a(t) * x^2 + b(t) * x
dfdx(x, t) = 2 * a(t) * x + b(t)

u(x,y,t) = -π * A * sin.(π * f(x, t)) .* cos.(π * y)
v(x,y,t) = π * A * cos.(π * f(x, t)) .* sin.(π * y) .* dfdx(x, t)
````

### Solve the Forward and Backward IVPs

````@example ftle_continuous
t0 = 0.0
T = 12
xf, yf = displacement_field(u,v,x0,y0,(t0,t0+T))
xb, yb = displacement_field(u,v,x0,y0,(t0,t0-T))
````

### Compute the Forward and Backward FTLEs

````@example ftle_continuous
fFTLE = similar(x0)
bFTLE = similar(x0)
compute_FTLE!(fFTLE,xf,yf,dx,dx,T)
compute_FTLE!(bFTLE,xb,yb,dx,dx,T)
````

### Plot the FTLEs on Top of Each Other

````@example ftle_continuous
plot(fFTLE, ftle_cache,fill=false, title="FTLE, t = $t0", xlabel="x", ylabel="y", colorbar=false, levels = 30, c=:inferno)
plot!(bFTLE, ftle_cache, fill=false, colorbar=false, levels = 30, c=:viridis)
````

The code here creates a gif

    fFTLE = similar(x0)
    bFTLE = similar(x0)

    @gif for t0 in 0.0:1.0:10.0
        print(t0)

        xf, yf = displacement_field(u,v,x0,y0,(t0,t0+T))
        xb, yb = displacement_field(u,v,x0,y0,(t0,t0-T))

        compute_FTLE!(fFTLE,xf,yf,dx,dx,T)
        compute_FTLE!(bFTLE,xb,yb,dx,dx,T)

        plot(fFTLE, ftle_cache,fill=false, title="FTLE, t = $t0", xlabel="x", ylabel="y", colorbar=false, levels = 30, c=:inferno)
        plot!(bFTLE, ftle_cache, fill=false, colorbar=false, levels = 30, c=:viridis)

    end every 1 fps = 2

## Example 2 - Issues with the Sliding-Window Approach
The sliding-window approach attempts to detect Langragian coherent structures (LCS) by computing the FTLE fields over windows of the form [t0, t0 + T] with varying t0 values. However, this approach does not obey Lagrangian invariance because the LCS's at different t0 values do not evolve into each other (Haller 2015).

This example illustrates the point above.

### First, generate initial conditions.

````@example ftle_continuous
X_MIN, X_MAX = -6.0, 50.0
Y_MIN, Y_MAX = -2, 2
dx = 0.08

ftlegrid = PhysicalGrid((X_MIN,X_MAX),(Y_MIN,Y_MAX),dx,optimize=false);
ftle_cache = SurfaceScalarCache(ftlegrid)
x0, y0 = x_grid(ftle_cache), y_grid(ftle_cache)
````

### Define the Vector Field

````@example ftle_continuous
u2(x,y) = 1 + tanh.(x).^2
v2(x,y) = -2 * tanh.(x) ./ cosh.(x).^2 .* y
````

### Solve the IVP and Compute FTLEs with the Sliding-Window Approach
As shown in the figure,although the flow is unsteady, the actual LCS is a material line that moves with the flow to the right. Nevertheless, as shown in the gif below, when t0 varies from 0.0 to 10.0, the forward FTLE fields reveals a repelling LCS fixed on the y-axis.

![Lack of Lagrangian Invariance](https://github.com/qiyuanbillwu/ILMPostProcessing.jl/raw/master/Lagrangian_Invariance.png)

````@example ftle_continuous
fFTLE = similar(x0)
T = 5

@gif for t0 in 0.0:0.5:10.0
    print(t0)

    xf, yf = displacement_field(u2,v2,x0,y0,(t0,t0+T))
    compute_FTLE!(fFTLE,xf,yf,dx,dx,T)

    plot(fFTLE, ftle_cache,fill=false, title="FTLE, t = $t0", xlabel="x", ylabel="y", colorbar=false, levels = 30, c=:inferno)

end every 1 fps = 2
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

