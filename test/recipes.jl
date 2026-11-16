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
        s = stability(pd; devs=(:adev, :hdev, :mdev), m=Octave, ci=false)
        @test only(RecipesBase.apply_recipe(Dict{Symbol,Any}(), s)).args == (s.results,)
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), s.results)
        @test length(rd) == 3
        # Decade ticks are set once for the axes, shared by every series.
        @test allequal(rdi.plotattributes[:yticks] for rdi in rd)
        @test all(rdi.plotattributes[:minorgrid] === true && !haskey(rdi.plotattributes, :subplot) for rdi in rd)
    end

    @testset "vector-of-results overlay" begin
        r0 = adev(pd, [1, 2, 4, 8]; ci=false)
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), [r, r0])
        @test length(rd) == 2
        # σ_y and the deviations in seconds get axes of their own.
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), [r, tdev(pd, [1, 2, 4, 8]), r0])
        pa = [rdi.plotattributes for rdi in rd]
        @test pa[1][:layout] == (1, 2) && [p[:subplot] for p in pa] == [1, 2, 1]
        @test pa[1][:left_margin] == pa[1][:bottom_margin] == (5, :mm)
        @test [p[:ylabel] for p in pa] == ["Deviation", "Deviation (s)", "Deviation"]
        @test pa[1][:yticks] == pa[3][:yticks]
        # The clocks of a hat take their names.
        h = ncornered_hat(ClockComparisons(randn(512, 2), 1.0; names=["a", "b"]), [1, 2, 4])
        @test [s.plotattributes[:label] for s in RecipesBase.apply_recipe(Dict{Symbol,Any}(), h)] == ["a", "b", "reference"]
    end

    @testset "dynamic deviation: 3-D waterfall" begin
        # The SP1065-style dynamic-deviation presentation: one σ(τ) curve per
        # window time, stacked along the time axis as 3-D lines. An
        # argument-less parent entry may also appear depending on the
        # RecipesBase version; the real curves are the 3-argument path3d ones.
        curves(rd) = [sd for sd in rd if length(sd.args) == 3 &&
                      get(sd.plotattributes, :seriestype, nothing) === :path3d]

        big = PhaseData(cumsum(randn(2048)) .* 1e-9, 1.0)
        res = dadev(big, [1, 2, 4]; window=128, step=512)   # 4 windows ⇒ labeled
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), res)   # default rendering
        cs = curves(rd)
        @test length(cs) == length(res.t) == 4        # one line per window
        attrs = [sd.plotattributes for sd in cs]
        @test all(pa -> pa[:seriestype] === :path3d && !haskey(pa, :yscale), attrs)
        # 3-D axes get no log scale, so τ and σ are plotted as log10; each
        # curve sits at its window time.
        @test all(cs[i].args[1] ≈ log10.(res.tau) && all(cs[i].args[2] .== res.t[i]) &&
                  cs[i].args[3] ≈ log10.(res.dev[i, :]) for i in eachindex(cs))
        @test [pa[:label] for pa in attrs] == ["t = $(round(t; sigdigits=4)) s" for t in res.t]
        # Decade ticks are integers labeled 1e<k> (TeX-safe in the docs build).
        @test all(all(t -> t == round(t), pa[k][1]) && all(startswith.(pa[k][2], "1e"))
                  for pa in attrs, k in (:xticks, :zticks))

        # Many-window maps suppress the per-curve legend.
        wide = dadev(big, [1, 2]; window=128, step=64)
        rdw = curves(RecipesBase.apply_recipe(Dict{Symbol,Any}(), wide))
        @test length(rdw) == length(wide.t) > 8
        @test all(sd -> sd.plotattributes[:label] == "", rdw)
    end

    @testset "spectral results" begin
        big = PhaseData(cumsum(randn(4096)) .* 1e-9, 1.0)
        # S_y / S_x render log-log with the DC bin (f = 0) dropped.
        for res in (psd_y(big; nperseg=512), psd_x(big; nperseg=512))
            rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), res)
            @test length(rd) == 1
            @test rd[1].plotattributes[:yscale] === :log10
            x, _ = rd[1].args
            @test minimum(x) > 0                 # no f = 0 on the log axis
            @test length(x) == length(res.freq) - 1
        end
        # ℒ(f) renders with a linear (dB) ordinate: decade ticks and minor
        # grid on the frequency axis only.
        rl = RecipesBase.apply_recipe(Dict{Symbol,Any}(), ssb_phase_noise(big; f_carrier=1e7, nperseg=512))
        @test length(rl) == 1
        @test !haskey(rl[1].plotattributes, :yscale)
        @test haskey(rl[1].plotattributes, :xticks)
        @test !haskey(rl[1].plotattributes, :yticks)
        @test rl[1].plotattributes[:xminorgrid] === true
        @test !haskey(rl[1].plotattributes, :minorgrid)
    end

    @testset "models, holdover, steering, time scales and filters" begin
        p = noise_gen(PhaseData, 4096, 1.0; h=Dict(0 => 2e-22, -2 => 1e-28), rng=Xoshiro(1))
        apply(x, attrs=Dict{Symbol,Any}()) = RecipesBase.apply_recipe(attrs, x)
        labels(rd) = [get(s.plotattributes, :label, "") for s in rd]
        # The panels a time-series recipe hands on, expanded to their lines.
        panels(x) = (rd = only(apply(x)); apply(only(rd.args), rd.plotattributes))

        fit = powerlaw_fit(adev(p); alphas=[0, -2])
        @test labels(apply(fit)) == ["measured", "fit", "WFM", "RWFM"]
        @test apply(fit)[3].args[2] ≈ powerlaw_sigma(Dict(0 => fit.h[1]), fit.x; tau0=1.0)
        sfit = powerlaw_fit(psd_y(p; nperseg=512); alphas=[2, 0, -2])
        @test sum(s.args[2] for s in apply(sfit)[3:end]) ≈ sfit.fit             # the components add up to the fit

        h = holdover(p; train=256, tau=[16, 64, 256], limit=1e-9)
        @test labels(apply(h)) == ["RMS", "95 % of trials within", "worst", "limit"]   # no Kalman σ for a fit
        @test last(apply(h)).plotattributes[:seriestype] === :hline
        hp = holdover(Dict(:WHFM => 2e-22, :RWFM => 1e-28); tau0=1.0, train=256, tau=[16, 64])
        @test labels(apply(hp)) == ["σ", "95 % within", "WFM", "RWFM"]

        m = clock_model(fit)
        g = steering_gains(300.0, 1.0)
        @test apply(steer(m, g; tau=[1.0, 10.0]))[1].args[1] == [1.0, 10.0]
        rs = panels(steer(p, m, g))
        @test [s.plotattributes[:subplot] for s in rs] == [1, 2] && rs[1].args[1][2] == 1.0

        kf = panels(kalman_filter(p, m))
        @test [s.plotattributes[:ylabel] for s in kf] == ["Phase (s)", "Frequency"]
        @test all(s -> length(s.plotattributes[:ribbon]) == 4096, kf)

        ms = [ClockModel(; tau0=1.0, q1=q) for q in (1e-24, 4e-24)]
        cc, _ = simulate_ensemble(ms, 512; reference=0, rng=Xoshiro(2))
        at1 = panels(timescale(cc; algorithm=:at1))
        @test [s.plotattributes[:subplot] for s in at1] == [1, 1, 2, 2]
        @test labels(at1) == ["clock 1", "clock 2", "", ""]                 # one legend
        @test length(panels(timescale(cc; algorithm=:jones_tryon, models=ms))) == 2      # no weights
    end
end
