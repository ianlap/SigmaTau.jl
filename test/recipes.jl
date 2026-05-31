# test/recipes.jl — RecipesBase extension smoke tests.
# Loading RecipesBase alongside SigmaTau triggers SigmaTauRecipesBaseExt.

using Test
using Random
using SigmaTau
using RecipesBase

@testset "Plot recipes" begin
    Random.seed!(99)
    pd = PhaseData(cumsum(randn(512)) .* 1e-9, 1.0)
    r  = adev(pd, [1, 2, 4, 8]; calc_ci=true)

    @testset "single result: error bars vs CI band" begin
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), r)
        @test length(rd) == 1
        @test haskey(rd[1].plotattributes, :yerror)
        @test !haskey(rd[1].plotattributes, :ribbon)

        rd_band = RecipesBase.apply_recipe(Dict{Symbol,Any}(:ci_band => true), r)
        @test haskey(rd_band[1].plotattributes, :ribbon)
        @test !haskey(rd_band[1].plotattributes, :yerror)
        @test !haskey(rd_band[1].plotattributes, :ci_band)   # consumed, not leaked
    end

    @testset "no CI ⇒ neither error bars nor ribbon" begin
        r0 = adev(pd, [1, 2, 4, 8]; calc_ci=false)
        rd0 = RecipesBase.apply_recipe(Dict{Symbol,Any}(), r0)
        @test !haskey(rd0[1].plotattributes, :yerror)
        @test !haskey(rd0[1].plotattributes, :ribbon)
    end

    @testset "suite overlay: one series per result" begin
        s = stability(pd; devs=(:adev, :hdev, :mdev), taus=Octave, calc_ci=false)
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), s)
        @test length(rd) == 3
    end

    @testset "vector-of-results overlay" begin
        r0 = adev(pd, [1, 2, 4, 8]; calc_ci=false)
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), [r, r0])
        @test length(rd) == 2
    end

    @testset "spectral results" begin
        big = PhaseData(cumsum(randn(4096)) .* 1e-9, 1.0)
        # S_y / S_x render log-log with the DC bin (f = 0) dropped.
        for res in (Sy(big; nperseg=512), Sx(big; nperseg=512))
            rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), res)
            @test length(rd) == 1
            @test rd[1].plotattributes[:yscale] === :log10
            x, _ = rd[1].args
            @test minimum(x) > 0                 # no f = 0 on the log axis
            @test length(x) == length(res.freq) - 1
        end
        # ℒ(f) renders with a linear (dB) ordinate.
        rl = RecipesBase.apply_recipe(Dict{Symbol,Any}(), L(big; f_carrier=1e7, nperseg=512))
        @test length(rl) == 1
        @test !haskey(rl[1].plotattributes, :yscale)
    end
end
