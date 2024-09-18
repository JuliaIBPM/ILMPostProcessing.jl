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

function euler_backward(initial_conditions, u, v, t0, t_start, dt, T)
    # this is not backward euler, it is forward euler method going back in time
    # t0 is the initial time 
    # t_start is where the starting time of the solution from viscous flow
    # dt is the interval between consecutive u and v fields
    # T is the integration time

    start_idx = Int(round((t0 - t_start) / dt + 1))
    iters = Int(round(T / dt - 1))

    z = initial_conditions

    for i = start_idx:-1:start_idx - iters
        x = z[:, 1]
        y = z[:, 2]
        z = z - dt * [u[i].(x, y) v[i].(x, y)]
    end

    return z

end

function compute_FTLE(FTLE, nx, ny, T, final_positions, dx, dy)

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