using ILMPostProcessing
using ViscousFlow
using LinearAlgebra

@testset "LAVD" begin

    # this test setup up an empty problem to compute LAVD fields on zero velocity fields
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

    # generate velocity fields and vorticity fields
    t_start = 0.0
    t_end = 1.0
    dt = 0.1
    tr = t_start:dt:t_end

    velxy = velocity_xy(sol,sys,tr)
    velseq = VectorFieldSequence(tr,velxy)
    vortxy = vorticity_xy(sol,sys,tr) 
    vortseq = ScalarFieldSequence(tr,vortxy)

    u, v = velxy[end]
    w = vortxy[end]

    @test norm(u) < 1e-8
    @test norm(v) < 1e-8
    @test norm(w) < 1e-8

    # generate initial conditions
    X_MIN = -0.5
    X_MAX = 0.5
    Y_MIN = -0.5
    Y_MAX = 0.5
    dx = 0.1
    lavdgrid = PhysicalGrid((X_MIN,X_MAX),(Y_MIN,Y_MAX),dx)
    lavd_cache = SurfaceScalarCache(lavdgrid)
    x0, y0 = x_grid(lavd_cache), y_grid(lavd_cache)

    # solve the forward and backward IVP using Euler and Adams-Bashforth
    # start at t = 0.5 and integrate to t= 1.0 in forward time and to t = 0.0 in backward time
    T = 0.5
    t0 = 0.5
    t1 = t0 + T

    # compute LAVD and IVD fields
    traj = compute_trajectory(velseq, (x0, y0), (t0, t1))
    LAVD = similar(x0)
    compute_LAVD!(LAVD, traj, X_MIN, Y_MIN, X_MAX, Y_MAX, vortseq)
    IVD = similar(x0)
    compute_IVD!(IVD, traj, X_MIN, Y_MIN, X_MAX, Y_MAX, vortseq)
    
    @test norm(LAVD) < 1e-8
    @test norm(IVD) < 1e-8

end