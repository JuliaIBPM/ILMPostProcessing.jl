using RecipesBase
using ColorTypes
#import PlotUtils: cgrad, palette, color_list

@recipe function f(traj::Trajectories; idxs=1:traj.np, trajname = "trajectory", startmarker=false)

  xguide --> "x"
  yguide --> "y"
  aspect_ratio := 1
  size --> (700,400)
  if !isa(idxs,AbstractVector{<:Int})
    error("idxs must be an AbstractVector of integers")
  end

  for jt in idxs
      x, y = traj[jt]
      @series begin
        linewidth --> 2
        label --> "$trajname $jt"
        x, y
      end
      if startmarker
        @series begin
          markersize --> 4
          markershape --> :circle
          label --> ""
          seriestype = :scatter
          [x[1]],[y[1]]
        end
      end
  end
end

@recipe function f(traj::Trajectories,field::T; idxs=1:traj.np, deriv=0, fieldname = "field") where T<:Union{Function,AbstractInterpolation,ScalarFieldSequence}

  if !isa(idxs,AbstractVector{<:Int})
    error("idxs must be an AbstractVector of integers")
  end
  for jt in idxs
    
      straj = field_along_trajectory(field,traj,jt,deriv=deriv)

      @series begin
        label --> "$fieldname on trajectory $jt"
        traj.t, straj
      end
    
  end
end


