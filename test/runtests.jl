using Test

@testset "SigmaTau" begin
    @testset "types"    begin include("types/runtests.jl") end
    @testset "stab"     begin include("stab/runtests.jl")  end
    @testset "io"       begin include("io/runtests.jl")    end
    @testset "umbrella" begin include("umbrella_smoke.jl") end
    @testset "recipes"  begin include("recipes.jl")        end
    @testset "tables"   begin include("tables.jl")         end
end
