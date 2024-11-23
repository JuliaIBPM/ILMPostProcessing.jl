"""
    compute_LAVD!(LAVD, traj, X_MIN, Y_MIN, X_MAX, Y_MAX, vortseq)

Compute the Lagrangian-averaged vorticity deviation `LAVD` field given the trajectories `traj` and the vorticity field sequence `vortseq`. 
"""
function compute_LAVD!(LAVD, traj, X_MIN, Y_MIN, X_MAX, Y_MAX, vortseq)
    t1 = traj.t[1]
    x1 = traj.xhistory[1]
    y1 = traj.yhistory[1]
    w3_1 = _instantaneous_scalar_field_in_series(vortseq, t1)  
    fn1 = abs.(w3_1.(x1, y1) - mean(w3_1.(x1, y1)))
    
    for i in 2:length(traj.t) 
        t2 = traj.t[i]
        x2 = traj.xhistory[i]
        y2 = traj.yhistory[i]
        w3_2 = _instantaneous_scalar_field_in_series(vortseq, t2)
        fn2 = abs.(w3_2.(x2, y2) - mean(w3_2.(x2, y2)))
        
        ds = sqrt.((x2-x1).^2+(y2-y1).^2)
        LAVD .= LAVD + 0.5 * (fn1 + fn2) .* ds

        t1 = t2
        x1 = x2
        y1 = y2
        fn1 = fn2
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
    IVD .= abs.(w3.(x0, y0) - mean(w3.(x0, y0)));
end

