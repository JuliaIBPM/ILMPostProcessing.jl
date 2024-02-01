@testset "POD" begin
    X = [rand(1000) for i=1:100]
    PODModes(X)
end