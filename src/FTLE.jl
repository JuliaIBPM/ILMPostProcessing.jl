# this function could be included in ViscousFlow.jl

function make_interp_fields!(u, v, t_start, t_end, dt, velocity_fn, sol, sys, grid)
    time = t_start:dt:t_end

    for t in time
        # Call the passed velocity function
        vel = velocity_fn(sol, sys, t)
        
        # Assuming vel contains u and v components
        push!(u, interpolatable_field(vel.u, grid))
        push!(v, interpolatable_field(vel.v, grid))
    end
end

function gen_init_conds(X_MIN, X_MAX, Y_MIN, Y_MAX, nx, ny)
    x0 = range(X_MIN, X_MAX, length=nx)
    y0 = range(Y_MIN, Y_MAX, length=ny)
    dx = x0[2] - x0[1]
    dy = y0[2] - y0[1]

    # Define the initial conditions as a standard matrix
    initial_conditions_matrix = [ [x0[i], y0[j]] for i in 1:nx, j in 1:ny]
    initial_conditions_matrix = initial_conditions_matrix'

    # Flatten the initial conditions matrix into a 1D array
    initial_conditions = vcat(initial_conditions_matrix...)
    
    return initial_conditions, dx, dy
end

function euler_forward(initial_conditions, u, v, t0, t_start, dt, T)
    # t0 is the initial time 
    # t_start is where the starting time of the solution from viscous flow
    # dt is the interval between consecutive u and v fields
    # T is the integration time

    w = initial_conditions

    if T == 0.0
        return w
    end

    start_idx = Int(round((t0 - t_start) / dt + 1))
    iters = Int(round(T / dt - 1))

    for i = start_idx:start_idx + iters
        x = w[:, 1]
        y = w[:, 2]
        w = w + dt * [u[i].(x, y) v[i].(x, y)]
    end

    return w

end

function euler_backward(initial_conditions, u, v, t0, t_start, dt, T)
    # this is not backward euler, it is forward euler method going back in time
    # t0 is the initial time 
    # t_start is where the starting time of the solution from viscous flow
    # dt is the interval between consecutive u and v fields
    # T is the integration time

    z = initial_conditions

    if T == 0.0
        return z
    end
    
    start_idx = Int(round((t0 - t_start) / dt + 1))
    iters = Int(round(T / dt - 1))

    for i = start_idx:-1:start_idx - iters
        x = z[:, 1]
        y = z[:, 2]
        z = z - dt * [u[i].(x, y) v[i].(x, y)]
    end

    return z

end

function adams_bashforth_2_forward(initial_conditions, u, v, t0, t_start, dt, T)
    # t0 is the initial time 
    # t_start is where the starting time of the solution from viscous flow
    # dt is the interval between consecutive u and v fields
    # T is the integration time

    if T == 0
        return initial_conditions
    end
    
    start_idx = Int(round((t0 - t_start) / dt + 1))
    iters = Int(round(T / dt - 1))

    # first generate the 2nd initial condition w2 by using forward euler once 
    w1 = initial_conditions
    x1 = w1[:, 1]
    y1 = w1[:, 2]
    w2 = w1 - dt * [u[start_idx].(x1, y1) v[start_idx].(x1, y1)]

    x2 = w2[:, 1]
    y2 = w2[:, 2]

    w3 = initial_conditions

    for i = 1:iters
        w3 = w2 + dt * (3/2*[u[start_idx+i].(x2, y2) v[start_idx+i].(x2, y2)] - 1/2*[u[start_idx+i-1].(x1, y1) v[start_idx+i-1].(x1, y1)])
        x1 = x2 
        y1 = y2
        x2 = w3[:,1]
        y2 = w3[:,2]
        w2 = w3
    end
    
    return w3
end

function adams_bashforth_2_backward(initial_conditions, u, v, t0, t_start, dt, T)
    # t0 is the initial time 
    # t_start is where the starting time of the solution from viscous flow
    # dt is the interval between consecutive u and v fields
    # T is the integration time

    if T == 0
        return initial_conditions
    end
    
    start_idx = Int(round((t0 - t_start) / dt + 1))
    iters = Int(round(T / dt - 1))

    # first generate the 2nd initial condition w2 by using forward euler once 
    w1 = initial_conditions
    x1 = w1[:, 1]
    y1 = w1[:, 2]
    w2 = w1 + dt * [u[start_idx].(x1, y1) v[start_idx].(x1, y1)]

    x2 = w2[:, 1]
    y2 = w2[:, 2]

    w3 = initial_conditions

    for i = -1:-1:-iters
        w3 = w2 - dt * (3/2*[u[start_idx+i].(x2, y2) v[start_idx+i].(x2, y2)] - 1/2*[u[start_idx+i+1].(x1, y1) v[start_idx+i+1].(x1, y1)])
        x1 = x2 
        y1 = y2
        x2 = w3[:,1]
        y2 = w3[:,2]
        w2 = w3
    end
    
    return w3
end

function compute_FTLE!(FTLE, nx, ny, T, final_positions, dx, dy)

    # Reshape the 2D array into the original 3D matrix format (nx, ny, 2)
    final_matrix = reshape(final_positions, ny, nx, 2)
    final_x = final_matrix[:,:,1];
    final_y = final_matrix[:,:,2];
    
    # Shifted arrays for vector operations
    final_x_i_minus = final_x[2:end-1, 1:end-2]
    final_x_i_plus  = final_x[2:end-1, 3:end]
    final_x_j_minus = final_x[1:end-2, 2:end-1]
    final_x_j_plus  = final_x[3:end, 2:end-1]

    final_y_i_minus = final_y[2:end-1, 1:end-2]
    final_y_i_plus  = final_y[2:end-1, 3:end]
    final_y_j_minus = final_y[1:end-2, 2:end-1]
    final_y_j_plus  = final_y[3:end, 2:end-1]

    # Compute the elements of the deformation gradient tensor A
    a11 = (final_x_i_plus - final_x_i_minus) / 2 / dx
    a12 = (final_x_j_plus - final_x_j_minus) / 2 / dy
    a21 = (final_y_i_plus - final_y_i_minus) / 2 / dx
    a22 = (final_y_j_plus - final_y_j_minus) / 2 / dy

    # Compute delta components
    delta_11 = a11.^2 .+ a21.^2
    delta_12 = a11 .* a12 .+ a21 .* a22
    delta_22 = a12.^2 .+ a22.^2

    # Eigenvalues of the delta matrix using characteristic equation
    trace = delta_11 .+ delta_22
    determinant = delta_11 .* delta_22 .- delta_12.^2
    discriminant = sqrt.(trace.^2 .- 4 .* determinant)
    lambda1 = (trace .+ discriminant) ./ 2

    # Compute FTLE (same slicing approach to match the dimensions)
    FTLE .= 1 / (2 * abs(T)) .* log.(lambda1)
end