using Test
using Random
using SigmaTau

@testset "SigmaTau IO" begin
    include("detrend.jl")
    include("fillgaps.jl")
    include("read.jl")
    include("results.jl")
    include("save.jl")
    include("outliers.jl")
end
