using Test
using Random
using SigmaTau

@testset "SigmaTau IO" begin
    include("detrend.jl")
    include("fillgaps.jl")
    include("read.jl")
end
