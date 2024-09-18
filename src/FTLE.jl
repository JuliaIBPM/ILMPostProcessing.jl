function euler_forward(initial_conditions, u, v, t0, t_start, dt, T)
    # t0 is the initial time 
    # t_start is where the starting time of the solution from viscous flow
    # dt is the interval between consecutive u and v fields
    # T is the integration time

    start_idx = Int(round((t0 - t_start) / dt + 1))
    iters = Int(round(T / dt - 1))

    w = initial_conditions

    for i = start_idx:start_idx + iters
        x = w[:, 1]
        y = w[:, 2]
        w = w + dt * [u[i].(x, y) v[i].(x, y)]
    end

    return w

end