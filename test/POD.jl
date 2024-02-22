using ILMPostProcessing
using ViscousFlow
using LinearAlgebra

@testset "POD" begin
    # setup grid
    xlim = (-1.0, 1.0)
    ylim = (-1.0, 1.0)
    g = PhysicalGrid(xlim, ylim, 0.01)
    # setup velocity field cache
    cache = SurfaceScalarCache(g)
    vel = zeros_gridgrad(cache)

    # create velocity snapshot data
    vsnap = [zeros_gridgrad(cache) for i=1:10]
    # create random velocity fields
    for v in vsnap
        v .= rand(size(v)...) # ... means splat 
    end

    # extract POD modes
    modes = PODModes(vsnap)
    # compare POD-reconstructed field with original vel field
    @test norm(modes.fieldReconst-vsnap[end])/norm(vsnap[end]) < 0.005
end