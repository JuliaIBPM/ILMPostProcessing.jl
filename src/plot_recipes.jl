using RecipesBase
using ColorTypes
#import PlotUtils: cgrad, palette, color_list

@recipe function f(traj::Trajectories)

  xguide --> "x"
  yguide --> "y"
  aspect_ratio := 1
  size --> (700,400)

  for jt in 1:traj.np
      @series begin
        linewidth --> 2
        label --> "particle $jt"
        traj[jt]
      end
  end
end

@recipe function f(traj::Trajectories,sys::ILMSystem)

  xguide --> "x"
  yguide --> "y"
  aspect_ratio := 1
  size --> (700,400)

  @series begin
    traj
  end
  @series begin
    sys.base_cache.bl
  end
  
end


@userplot FieldTrajectory

@recipe function f(h::FieldTrajectory;fieldlabel="Field",deriv=0)
  if length(h.args) != 3
      error("`fieldtrajectory` should be given three arguments.  Got: $(typeof(h.args))")
  end
  traj, field, sys = h.args

  xtraj, ytraj = traj[1]

  layout := (2,1)
  size --> (600,600)
  xlims -> (-1,2)

  @series begin
    subplot := 1
    aspect_ratio := 1
    ylims --> (-0.5,0.5)
    title := "Particle trajectory"
    traj
  end

  if length(sys) > 0
    body = surfaces(sys)[1]
    @series begin
      subplot := 1
      ylims --> (-0.5,0.5)
      xguide --> "x"
      yguide --> "y"
      aspect_ratio := 1
      surfaces(sys)
    end
    yb = -20
    yt = 20
    Xr = [minimum(body.x), maximum(body.x), maximum(body.x), minimum(body.x), minimum(body.x)]
    Yr = [yb,yb,yt,yt,yb]
    @series begin
      subplot := 2
      fillrange := 0
      fillalpha --> 0.2
      linealpha --> 0.2
      fillcolor := :lightgray
      linecolor := :lightgray
      label := "Body location"
      Xr,Yr
    end
  end

  
  if typeof(field) <: VectorGridData
    utraj,vtraj = field_along_trajectory(field,sys,traj,deriv=deriv)
    minuv = min(minimum(utraj),minimum(vtraj)) - 0.5
    maxuv = max(maximum(utraj),maximum(vtraj)) + 0.5

    @series begin
      subplot := 2
      xguide := "x"
      label := "x %$fieldlabel"
      xtraj, utraj
    end

    @series begin
      subplot := 2
      xguide := "x"
      label := "y %$fieldlabel"
      title := "$fieldlabel components along trajectory"
      legend := true
      ylims := (minuv,maxuv)
      xtraj, vtraj
    end
  elseif typeof(field) <: ScalarGridData
    straj = field_along_trajectory(field,sys,traj,deriv=deriv)

    @series begin
      subplot := 2
      xguide := "x"
      label := "$fieldlabel"
      title := "$fieldlabel along trajectory"
      ylims := (minimum(straj)-0.5,maximum(straj)+0.5)
      legend := true
      xtraj, straj
    end
  end

end


