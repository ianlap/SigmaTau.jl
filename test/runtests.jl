using Test

@testset "SigmaTau" begin
    @testset "types"    begin include("types/runtests.jl") end
    @testset "stab"     begin include("stab/runtests.jl")  end
    @testset "est"      begin include("est/runtests.jl")   end
    @testset "umbrella" begin include("umbrella_smoke.jl") end
end
