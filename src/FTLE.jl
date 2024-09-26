"""
    make_interp_fields!(u, v, t_start, t_end, dt, velocity_fn, sol, sys, grid)

Generate an array of interpolatable velocity fields u and v using the solution of a ViscousFlow problem. 

Each element of `u` and `v` is essentially a "Vectorfield!" used when solving an ODEProblem. Note: This function could be included in ViscousFlow.jl. 

# Arguments
- `u`: x-velocity fields
- `v`: y-velocity fields
- `t_start`: start time of sol
- `t_end`: end time of sol 
- `dt`: step size between consecutive velocity fields
- `velocity_fn`: function to compute velocity from solution, should be ViscousFlow.velocity
- `sys`: viscous flow system
- `grid`: physical grid
"""
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

"""
    gen_init_conds(X_MIN, X_MAX, Y_MIN, Y_MAX, nx, ny)

Generate a list of initial points (x, y).

These initial conditions represent a collocated grid with `nx` points in (`X_MIN`, `X_MAX`) and `ny` points in (`Y_MIN`, `Y_MAX`). The points are then flattened to a 1D array. The initial conditions could be used to compute FTLE or to visaulize trajectories.
"""
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

"""
    euler_forward(initial_conditions, u, v, t0, t_start, dt, T)

Solve the initial value problem (IVP) using forward Euler method.

Integrate forward in time to compute forward FTLE fields.

# Arguments
- `initial_conditions`: generated with function gen_init_conds
- `u`: array of x-velocity fields
- `v`: array of y-velocity fields
- `t0`: initial time
- `t_start`: start time of `u` and `v`
- `dt`: step size between consecutive velocity fields
- `T`: length of integration time
"""
function euler_forward(initial_conditions, u, v, t0, t_start, dt, T)
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

"""
    euler_backward(initial_conditions, u, v, t0, t_start, dt, T)

Solve the initial value problem (IVP) using forward Euler method.

Integrate backward in time to compute backward FTLE fields. Note: not backward Euler method.

# Arguments
- `initial_conditions`: generated with function gen_init_conds
- `u`: array of x-velocity fields
- `v`: array of y-velocity fields
- `t0`: initial time
- `t_start`: start time of `u` and `v`
- `dt`: step size between consecutive velocity fields
- `T`: length of integration time
"""
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

"""
    adams_bashforth_2_forward(initial_conditions, u, v, t0, t_start, dt, T)

Solve the initial value problem (IVP) using 2-step Adams-Bashforth method.

Integrate forward in time to compute forward FTLE fields.

# Arguments
- `initial_conditions`: generated with function gen_init_conds
- `u`: array of x-velocity fields
- `v`: array of y-velocity fields
- `t0`: initial time
- `t_start`: start time of `u` and `v`
- `dt`: step size between consecutive velocity fields
- `T`: length of integration time
"""
function adams_bashforth_2_forward(initial_conditions, u, v, t0, t_start, dt, T)
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

"""
    adams_bashforth_2_backward(initial_conditions, u, v, t0, t_start, dt, T)

Solve the initial value problem (IVP) using 2-step Adams-Bashforth method.

Integrate backward in time to compute backward FTLE fields.

# Arguments
- `initial_conditions`: generated with function gen_init_conds
- `u`: array of x-velocity fields
- `v`: array of y-velocity fields
- `t0`: initial time
- `t_start`: start time of `u` and `v`
- `dt`: step size between consecutive velocity fields
- `T`: length of integration time
"""
function adams_bashforth_2_backward(initial_conditions, u, v, t0, t_start, dt, T)
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

"""
    compute_FTLE!(FTLE, nx, ny, T, final_positions, dx, dy)

Compute the `FTLE` field given the final positions of initial points on a collocated grid. 

The underlying math is detailed in: https://shaddenlab.berkeley.edu/uploads/LCS-tutorial/computation.html. For each grid point, first compute the gradient of the flow map using two point central differencing. Then, calculate the maximum eigenvalue of the 2 x 2 gradient matrix. Finally, compute the FTLE value using the eigenvalue.

# Arguments
- `FTLE`: an empty 2D array (i.e., `FTLE = zeros(Float64, ny - 2, nx - 2))`, `nx - 2` and `ny - 2` accounts for the boundary points in the central difference formula
- `nx`: number of grid points in x 
- `ny`: number of grid points in y 
- `T`: length of integration time
- `final_positions`: solutions of the IVP 
- `dx`: spacing of initial grids in x 
- `dy`: spacing of initial grids in y 
"""
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

    # Compute the components of delta matrix = A' * A
    a = a11.^2 .+ a21.^2
    b = a11 .* a12 .+ a21 .* a22
    c = a12.^2 .+ a22.^2

    # Eigenvalues of the delta matrix using characteristic equation
    lambda = (a .+ c .+ sqrt.((a .- c).^2 .+ 4 .* b.^2)) ./ 2

    # Compute FTLE (same slicing approach to match the dimensions)
    FTLE .= 1 / (2 * abs(T)) .* log.(lambda)
end


"""
    compute_FTLE!(FTLE::ScalarGridData, final_x::ScalarGridData, final_y::ScalarGridData, dx::Real, dy::Real, T::Real)

Compute the `FTLE` field given the final positions of initial points on a collocated grid. 

The underlying math is detailed in: https://shaddenlab.berkeley.edu/uploads/LCS-tutorial/computation.html. For each grid point, first compute the gradient of the flow map using two point central differencing. Then, calculate the maximum eigenvalue of the 2 x 2 gradient matrix. Finally, compute the FTLE value using the eigenvalue.

# Arguments
- `FTLE`: will hold the FTLE field, in `ScalarGridData` type
- `final_x`: deformed x positions, in `ScalarGridData` type
- `final_y`: deformed y positions, in `ScalarGridData` type 
- `dx`: spacing of initial grids in x 
- `dy`: spacing of initial grids in y 
- `T`: length of integration time
"""
function compute_FTLE!(FTLE::ScalarGridData, final_x::ScalarGridData, final_y::ScalarGridData, dx::Real, dy::Real, T::Real)

    nx, ny = size(final_x)
    # Shifted arrays for vector operations
    final_x_i_minus = view(final_x,2:nx-1, 1:ny-2)
    final_x_i_plus  = view(final_x,2:nx-1, 3:ny)
    final_x_j_minus = view(final_x,1:nx-2, 2:ny-1)
    final_x_j_plus  = view(final_x,3:nx, 2:ny-1)

    final_y_i_minus = view(final_y,2:nx-1, 1:ny-2)
    final_y_i_plus  = view(final_y,2:nx-1, 3:ny)
    final_y_j_minus = view(final_y,1:nx-2, 2:ny-1)
    final_y_j_plus  = view(final_y,3:nx, 2:ny-1)

    # Compute the elements of the deformation gradient tensor A
    a11 = (final_x_i_plus - final_x_i_minus) / 2 / dx
    a12 = (final_x_j_plus - final_x_j_minus) / 2 / dy
    a21 = (final_y_i_plus - final_y_i_minus) / 2 / dx
    a22 = (final_y_j_plus - final_y_j_minus) / 2 / dy

    # Compute the components of delta matrix = A' * A
    a = a11.^2 .+ a21.^2
    b = a11 .* a12 .+ a21 .* a22
    c = a12.^2 .+ a22.^2

    # Eigenvalues of the delta matrix using characteristic equation
    lambda = (a .+ c .+ sqrt.((a .- c).^2 .+ 4 .* b.^2)) ./ 2

    # Compute FTLE (same slicing approach to match the dimensions)
    FTLE_center = view(FTLE,2:nx-1,2:ny-1)
    FTLE_center .= 1 / (2 * abs(T)) .* log.(lambda)
end