using ILMPostProcessing
using ViscousFlow
using LinearAlgebra

@testset "FTLE" begin

    # this test setup up an empty problem to compute FTLE fields on zero velocity fields
    # setup and solve viscous flow problem
    my_params = Dict()
    my_params["Re"] = 200
    xlim = (-1.0,1.0)
    ylim = (-1.0,1.0)

    g = setup_grid(xlim,ylim,my_params)
    sys = viscousflow_system(g,phys_params=my_params)

    u0 = init_sol(sys)

    tspan = (0.0,1.0)
    integrator = init(u0,tspan,sys)
    step!(integrator,1.0)
    sol = integrator.sol

    # generate velocity fields
    #u = []
    #v = []
    t_start = 0.0
    t_end = 1.0
    dt = 0.1
    tr = t_start:dt:t_end

    velxy = velocity_xy(sol,sys,tr)
    vseq = VectorFieldSequence(tr,velxy)

    u, v = velxy[end]

    @test norm(u) < 1e-8
    @test norm(v) < 1e-8

    # generate initial conditions
    X_MIN = -0.5
    X_MAX = 0.5
    Y_MIN = -0.5
    Y_MAX = 0.5
    #nx, ny = 10, 10
    dx = 0.1
    ftlegrid = PhysicalGrid((X_MIN,X_MAX),(Y_MIN,Y_MAX),dx)
    ftle_cache = SurfaceScalarCache(ftlegrid)

    #initial_conditions, dx, dy = ILMPostProcessing.gen_init_conds(X_MIN, X_MAX, Y_MIN, Y_MAX, nx, ny)
    x0, y0 = x_grid(ftle_cache), y_grid(ftle_cache)

    # solve the forward and backward IVP using Euler and Adams-Bashforth
    # start at t = 0.5 and integrate to t= 1.0 in forward time and to t = 0.0 in backward time
    T = 0.5
    t0 = 0.5

    xf, yf = displacement_field(vseq,x0,y0,(t0,t0+T),alg=ILMPostProcessing.Euler())
    xb, yb = displacement_field(vseq,x0,y0,(t0,t0-T),alg=ILMPostProcessing.Euler())

    #a = ILMPostProcessing.euler_forward(initial_conditions, u, v, t0, t_start, dt, T)
    #b = ILMPostProcessing.euler_backward(initial_conditions, u, v, t0, t_start, dt, T)
    #c = ILMPostProcessing.adams_bashforth_2_forward(initial_conditions, u, v, t0, t_start, dt, T)
    #d = ILMPostProcessing.adams_bashforth_2_backward(initial_conditions, u, v, t0, t_start, dt, T)

    #@test norm(a - initial_conditions) < 1e-8
    #@test norm(b - initial_conditions) < 1e-8
    #@test norm(c - initial_conditions) < 1e-8
    #@test norm(d - initial_conditions) < 1e-8

    @test norm(xf - x0) < 1e-8
    @test norm(yf - y0) < 1e-8
    @test norm(xb - x0) < 1e-8
    @test norm(yf - y0) < 1e-8

    # compute FTLE fields, should all be zeros
    #ftle_a = ftle_b = ftle_c = ftle_d = zeros(Float64, ny - 2, nx - 2)
    #ILMPostProcessing.compute_FTLE!(ftle_a, nx, ny, T, a, dx, dy)
    #ILMPostProcessing.compute_FTLE!(ftle_b, nx, ny, T, b, dx, dy)
    #ILMPostProcessing.compute_FTLE!(ftle_c, nx, ny, T, c, dx, dy)
    #ILMPostProcessing.compute_FTLE!(ftle_d, nx, ny, T, d, dx, dy)
    fFTLE = similar(x0)
    compute_FTLE!(fFTLE,xf,yf,dx,dx,T)

    bFTLE = similar(x0)
    compute_FTLE!(bFTLE,xb,yb,dx,dx,T)

    @test norm(fFTLE) < 1e-8
    @test norm(bFTLE) < 1e-8
    #@test norm(ftle_c) < 1e-8
    #@test norm(ftle_d) < 1e-8

end