"""
    _integrate_w3(t, X_MIN, Y_MIN, X_MAX, Y_MAX, vortseq)

Calculate the mean vorticity in the specified domain. 
"""
function _integrate_w3(t, X_MIN, Y_MIN, X_MAX, Y_MAX, vortseq)
    w3 = _instantaneous_scalar_field_in_series(vortseq, t)
    integrand = x -> w3(x[1], x[2])  # x[1] = x, x[2] = y
    result, error_estimate = hcubature(integrand, [X_MIN, Y_MIN], [X_MAX, Y_MAX])  
    return result / ((X_MAX - X_MIN) * (Y_MAX - Y_MIN))
end

"""
    compute_LAVD!(LAVD, traj, X_MIN, Y_MIN, X_MAX, Y_MAX, vortseq)

Compute the Lagrangian-averaged vorticity deviation `LAVD` field given the trajectories `traj` and the vorticity field sequence `vortseq`. 
"""
function compute_LAVD!(LAVD, traj, X_MIN, Y_MIN, X_MAX, Y_MAX, vortseq)
    for i in length(traj.t) - 1
        t1 = traj.t[i]
        t2 = traj.t[i+1]
        x1 = traj.xhistory[i]
        x2 = traj.xhistory[i+1]
        y1 = traj.yhistory[i]
        y2 = traj.yhistory[i+1]
    
        w3_1 = _instantaneous_scalar_field_in_series(vortseq, t1)
        w3_2 = _instantaneous_scalar_field_in_series(vortseq, t2)
    
        fn1 = abs.(w3_1.(x1, y1) - _integrate_w3(t1, X_MIN, Y_MIN, X_MAX, Y_MAX, vortseq))
        fn2 = abs.(w3_2.(x2, y2) - _integrate_w3(t2, X_MIN, Y_MIN, X_MAX, Y_MAX, vortseq))
        
        ds = sqrt.((x2-x1).^2+(y2-y1).^2)
        LAVD .= LAVD + 0.5 * (fn1 + fn2) .* ds
    end
end

"""
    compute_IVD!(IVD, traj, X_MIN, Y_MIN, X_MAX, Y_MAX, vortseq)

Compute the instantaneous vorticity deviation `IVD` field given the trajectories `traj` and the vorticity field sequence `vortseq`. 
"""
function compute_IVD!(IVD, traj, X_MIN, Y_MIN, X_MAX, Y_MAX, vortseq)
    t0 = traj.t[1]
    x0 = traj.xhistory[1]
    y0 = traj.yhistory[1]
    w3 = _instantaneous_scalar_field_in_series(vortseq, t0)
    IVD .= abs.(w3.(x0, y0) - _integrate_w3(t0, X_MIN, Y_MIN, X_MAX, Y_MAX, vortseq));
end

