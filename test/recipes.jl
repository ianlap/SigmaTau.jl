# test/recipes.jl — RecipesBase extension smoke tests.
# Loading RecipesBase alongside SigmaTau triggers SigmaTauRecipesBaseExt.

using Test
using Random
using SigmaTau
using RecipesBase

@testset "Plot recipes" begin
    Random.seed!(99)
    pd = PhaseData(cumsum(randn(512)) .* 1e-9, 1.0)
    r  = adev(pd, [1, 2, 4, 8]; ci=true)

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

    @testset "marker and decade-tick defaults" begin
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), r)
        pa = rd[1].plotattributes
        @test pa[:markershape] === :circle
        @test pa[:markersize] == 3
        @test pa[:minorgrid] === true
        @test pa[:minorgridalpha] < pa[:gridalpha]   # minor fainter than major
        # Ticks sit at integer powers of 10 and span the data range.
        for ticks in (pa[:xticks], pa[:yticks])
            @test all(t -> t ≈ 10.0^round(log10(t)), ticks)
        end
        @test first(pa[:xticks]) <= minimum(r.tau) <= maximum(r.tau) <= last(pa[:xticks])
        @test first(pa[:yticks]) <= minimum(ci_lower(r)) <= maximum(ci_upper(r)) <= last(pa[:yticks])
        # User attributes override every default.
        rd2 = RecipesBase.apply_recipe(
            Dict{Symbol,Any}(:markershape => :square, :xticks => [1.0, 10.0]), r)
        @test rd2[1].plotattributes[:markershape] === :square
        @test rd2[1].plotattributes[:xticks] == [1.0, 10.0]
    end

    @testset "no CI ⇒ neither error bars nor ribbon" begin
        r0 = adev(pd, [1, 2, 4, 8]; ci=false)
        rd0 = RecipesBase.apply_recipe(Dict{Symbol,Any}(), r0)
        @test !haskey(rd0[1].plotattributes, :yerror)
        @test !haskey(rd0[1].plotattributes, :ribbon)
    end

    @testset "suite overlay: one series per result" begin
        s = stability(pd; devs=(:adev, :hdev, :mdev), taus=Octave, ci=false)
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), s)
        @test length(rd) == 3
        # Decade ticks are set once at the suite level, shared by every series.
        @test all(haskey(rdi.plotattributes, :xticks) for rdi in rd)
        @test all(rdi.plotattributes[:minorgrid] === true for rdi in rd)
    end

    @testset "vector-of-results overlay" begin
        r0 = adev(pd, [1, 2, 4, 8]; ci=false)
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), [r, r0])
        @test length(rd) == 2
    end

    @testset "dynamic deviation map: log10 heatmap" begin
        big = PhaseData(cumsum(randn(2048)) .* 1e-9, 1.0)
        res = dadev(big, [1, 2, 4, 100]; window=128)   # m=100 ⇒ NaN column
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), res)
        # All series are emitted via @series; an argument-less parent entry may
        # also appear depending on the RecipesBase version. Keep the data series.
        data = [sd for sd in rd if length(sd.args) == 3]
        @test length(data) == 1
        pa = data[1].plotattributes
        @test pa[:seriestype] === :heatmap
        @test pa[:yscale] === :log10                  # τ axis is log
        @test haskey(pa, :colorbar_title)             # log10 σ legend
        @test occursin("log10", string(pa[:colorbar_title]))
        # Decade ticks on the τ axis span the grid.
        @test all(t -> t ≈ 10.0^round(log10(t)), pa[:yticks])
        # Heatmap args: x = t, y = τ, z = log10(dev) transposed to τ × t.
        x, y, z = data[1].args
        @test x == res.t && y == res.tau
        @test size(z) == (length(res.tau), length(res.t))
        @test z[1, 1] ≈ log10(res.dev[1, 1])
        @test all(isnan, z[end, :])                   # unsupported m stays NaN
    end

    @testset "dynamic deviation map: 3-D waterfall mode" begin
        # The SP1065-style dynamic-deviation presentation: one σ(τ) curve per
        # window time, stacked along the time axis as 3-D lines. An
        # argument-less parent entry may also appear depending on the
        # RecipesBase version; the real curves are the 3-argument path3d ones.
        curves(rd) = [sd for sd in rd if length(sd.args) == 3 &&
                      get(sd.plotattributes, :seriestype, nothing) === :path3d]

        big = PhaseData(cumsum(randn(2048)) .* 1e-9, 1.0)
        res = dadev(big, [1, 2, 4]; window=128, step=512)   # 4 windows ⇒ labelled
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(:seriestype => :path3d), res)
        cs = curves(rd)
        @test length(cs) == length(res.t) == 4        # one line per window
        for (i, sd) in enumerate(cs)
            pa = sd.plotattributes
            @test pa[:seriestype] === :path3d
            @test !haskey(pa, :yscale)                # 3-D axes get no log scale…
            x, y, z = sd.args
            @test x ≈ log10.(res.tau)                 # …so τ is plotted as log10(τ)
            @test all(y .== res.t[i])                 # the curve sits at its t
            @test z ≈ log10.(res.dev[i, :])           # σ plotted as log10(σ)
            @test pa[:label] == "t = $(round(res.t[i]; sigdigits=4)) s"
            # Decade tick positions are integers labelled in plain-text 1e<k>
            # form (TeX-safe under the PGFPlotsX docs build).
            for key in (:xticks, :zticks)
                ticks, labels = pa[key]
                @test all(t -> t == round(t), ticks)
                @test all(startswith.(labels, "1e"))
            end
        end

        # Many-window maps suppress the per-curve legend.
        wide = dadev(big, [1, 2]; window=128, step=64)
        rdw = curves(RecipesBase.apply_recipe(Dict{Symbol,Any}(:seriestype => :path3d), wide))
        @test length(rdw) == length(wide.t) > 8
        @test all(sd -> sd.plotattributes[:label] == "", rdw)
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
        # ℒ(f) renders with a linear (dB) ordinate: decade ticks and minor
        # grid on the frequency axis only.
        rl = RecipesBase.apply_recipe(Dict{Symbol,Any}(), L(big; f_carrier=1e7, nperseg=512))
        @test length(rl) == 1
        @test !haskey(rl[1].plotattributes, :yscale)
        @test haskey(rl[1].plotattributes, :xticks)
        @test !haskey(rl[1].plotattributes, :yticks)
        @test rl[1].plotattributes[:xminorgrid] === true
        @test !haskey(rl[1].plotattributes, :minorgrid)
    end
end
