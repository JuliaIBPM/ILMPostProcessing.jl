using ILMPostProcessing
using ViscousFlow
using LinearAlgebra

@testset "Trajectories from computed flow fields" begin
    my_params = Dict()
    my_params["Re"] = 100
    xlim = (-2.0,2.0)
    ylim = (-2.0,2.0)
    my_params["grid Re"] = 10.0
    g = setup_grid(xlim,ylim,my_params)
    sys = viscousflow_system(g,phys_params=my_params)

    σ = 0.1
    x01, y01 = 1.0, 0.0
    x02, y02 = -1.0, 0.0
    A = 3
    twogauss = SpatialGaussian(σ,σ,x01,y01,A) + SpatialGaussian(σ,σ,x02,y02,A)
    u0 = init_sol(twogauss,sys)

    T = 1.0
    tspan = (0.0,T)
    integrator = init(u0,tspan,sys)

    step!(integrator,1.0)

    sol = integrator.sol

    velfcn = velocity_xy(integrator)

    pts = [ [0.5,1], [-0.5,-0.5], [0.25, -0.3]]

    traj = compute_trajectory(velfcn...,pts,(0,1),alg=Tsit5())

    @test length(traj.xhistory[end]) == length(traj.yhistory[end]) == length(pts)
   
    tsline = integrator.t
    sline = compute_streamline(velfcn...,pts,(0,12),tsline,alg=Tsit5())


    t_start = 0.0
    t_end = 1.0
    dt = timestep(u0,sys)
    tr = t_start:dt:t_end
    velxy = velocity_xy(sol,sys,tr)
    vseq = VectorFieldSequence(tr,velxy)

    dx = 0.1
    ftlegrid = PhysicalGrid((-2,2),(-2,2),dx)
    ftle_cache = SurfaceScalarCache(ftlegrid)
    x0, y0 = x_grid(ftle_cache), y_grid(ftle_cache)
    
    traj = compute_trajectory(vseq,(x0,y0),(t_start,t_end),alg=Euler())

    @test length(traj.xhistory[end]) == length(traj.yhistory[end]) == length(x0)


end


@testset "Trajectories from velocity functions" begin

    ufcn(x,y) = x
    vfcn(x,y) = -y

    pts = [ [0.5,1], [-0.5,-0.5], [0.25, -0.3]]

    traj = compute_trajectory(ufcn,vfcn,pts,(0,1),alg=Tsit5())

    @test typeof(traj.xhistory)<:Vector && typeof(traj.yhistory)<:Vector

    ufcn2(x,y,t) = 1
    vfcn2(x,y,t) = sin(t)

    traj = compute_trajectory(ufcn2,vfcn2,pts,(0,1),alg=Tsit5())

    @test typeof(traj.xhistory)<:Vector && typeof(traj.yhistory)<:Vector

    x1, y1 = traj[1]
    @test norm(pts[1][2] .+ 1.0 .- cos.(traj.t) .- y1) < 1e-7


    tsline = π
    sline = compute_streamline(ufcn2,vfcn2,pts,(0,10),tsline,alg=Tsit5())
    x1, y1 = sline[2]
    @test norm(y1 .- pts[2][2]) < 1e-7

    u(x,y,t) = 1.0
    v(x,y,t) = cos(2π*(x-t))
    y = [0.0,0.0]
    t = 0.0
    streak = compute_streakline(u,v,y,t)
    xs, ys = streak[1]
    @test xs[1] ≈ 3 && ys[1] ≈ 3 && xs[end] ≈ y[1] && ys[end] ≈ y[2]



end
