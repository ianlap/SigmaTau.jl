# test/stab/dynamic.jl — Dynamic (time-resolved) deviations: DADEV / DHDEV.
#
# Both estimators slide an analysis window across the record and evaluate the
# existing overlapping kernels (_adev_core / _hdev_core) per window, so the
# tests anchor against the static estimators (full-record window identity),
# the NaN support pattern, and the time localization that is the whole point
# of the dynamic construction.

@testset "Dynamic deviations (DADEV / DHDEV)" begin
    @testset "result geometry and window-center times" begin
        Random.seed!(20260610)
        N      = 2048
        window = 512
        step   = 256
        pd = PhaseData(cumsum(randn(N)) .* 1e-9, 1.0)
        ms = [1, 2, 4, 8]

        r = dadev(pd, ms; window=window, step=step)
        @test r isa DynamicStabilityResult
        @test r.deviation_type === :dadev
        @test r.window == window
        @test r.tau0 == 1.0
        @test r.tau == Float64.(ms)
        # Window starts 1, 1+step, … while start + window − 1 ≤ N.
        n_windows = length(1:step:(N - window + 1))
        @test size(r.dev) == (n_windows, length(ms))
        @test length(r.t) == n_windows
        # Window covering samples s … s+window−1 is centered at
        # (s − 1 + (window − 1)/2)·τ₀.
        @test r.t[1] ≈ (window - 1) / 2
        @test r.t[2] - r.t[1] ≈ step * 1.0
        @test all(isfinite, r.dev)

        # step default is window ÷ 2.
        r_def = dadev(pd, ms; window=window)
        @test r_def.t == r.t
        @test r_def.dev == r.dev
    end

    @testset "stationary record ⇒ near-constant dev along t" begin
        # White FM at a fixed level: every window sees the same process, so
        # the map should be flat along t up to estimator scatter. ADEV's
        # relative spread at (window, m) is ≈ 1/√(window − 2m); with
        # window = 1024 and m ≤ 4 that is ≈ 3 %, so a ±25 % band around the
        # per-τ median is generous.
        Random.seed!(42)
        N      = 8192
        window = 1024
        pd = PhaseData(_gen_powerlaw_phase(0.0, N; tau0=1.0), 1.0)
        ms = [1, 2, 4]

        for r in (dadev(pd, ms; window=window), dhdev(pd, ms; window=window))
            @test all(isfinite, r.dev)
            for j in eachindex(ms)
                col = r.dev[:, j]
                mid = sort(col)[cld(length(col), 2)]   # median
                @test all(c -> 0.75 * mid < c < 1.25 * mid, col)
            end
        end
    end

    @testset "noise step is localized at the right t" begin
        # First half white FM at level σ, second half at 10σ: windows fully
        # inside each half read the corresponding level, so the map shows the
        # step at t ≈ N/2 · τ₀. The 10× variance ratio leaves lots of margin
        # over the ≈3 % estimator scatter at window = 512, m = 1.
        Random.seed!(7)
        N      = 4096
        window = 512
        step   = 128
        y = randn(N - 1)
        y[(N ÷ 2):end] .*= 10.0
        pd = PhaseData(pushfirst!(cumsum(y), 0.0), 1.0)

        r = dadev(pd, [1]; window=window, step=step)
        col = r.dev[:, 1]
        t_step = (N ÷ 2) * 1.0                       # step time in seconds
        early = col[r.t .<= t_step - window / 2]     # windows fully pre-step
        late  = col[r.t .>= t_step + window / 2]     # windows fully post-step
        @test !isempty(early) && !isempty(late)
        # Levels are right (10× apart) and each plateau is internally flat.
        @test minimum(late) > 5 * maximum(early)
        @test maximum(early) < 1.5 * minimum(early)
        @test maximum(late)  < 1.5 * minimum(late)
        # The transition is bracketed: the first window touching the step is
        # the first to read above the early plateau.
        first_high = findfirst(c -> c > 2 * maximum(early), col)
        @test first_high !== nothing
        @test r.t[first_high] > t_step - window      # not before any overlap
        @test r.t[first_high] < t_step + window      # not after full overlap
    end

    @testset "full-record window reproduces the static estimators" begin
        Random.seed!(99)
        N  = 1024
        pd = PhaseData(cumsum(randn(N)) .* 1e-9, 1.0)
        ms = [1, 2, 4, 8, 16]

        # window = N admits exactly one window regardless of step ≥ 1, and
        # the single row is the plain overlapping estimate on the record.
        rd = dadev(pd, ms; window=N, step=N)
        @test size(rd.dev, 1) == 1
        @test vec(rd.dev) == adev(pd, ms; ci=false).dev

        rh = dhdev(pd, ms; window=N, step=2N)
        @test size(rh.dev, 1) == 1
        @test vec(rh.dev) == hdev(pd, ms; ci=false).dev
    end

    @testset "NaN pattern at unsupported m (per window, not per record)" begin
        Random.seed!(3)
        N      = 2048
        window = 64
        pd = PhaseData(cumsum(randn(N)) .* 1e-9, 1.0)

        # DADEV needs window − 2m ≥ 2: m = 31 works at window = 64, 32 not.
        rd = dadev(pd, [1, 31, 32, 100]; window=window)
        @test all(isfinite, rd.dev[:, 1])
        @test all(isfinite, rd.dev[:, 2])
        @test all(isnan,    rd.dev[:, 3])
        @test all(isnan,    rd.dev[:, 4])

        # DHDEV needs window − 3m ≥ 2: m = 20 works at window = 64, 21 not.
        rh = dhdev(pd, [1, 20, 21]; window=window)
        @test all(isfinite, rh.dev[:, 1])
        @test all(isfinite, rh.dev[:, 2])
        @test all(isnan,    rh.dev[:, 3])
    end

    @testset "TauMode grids clamp to the window, not the record" begin
        Random.seed!(11)
        N      = 4096
        window = 128
        pd = PhaseData(cumsum(randn(N)) .* 1e-9, 1.0)

        rd = dadev(pd, Octave; window=window)
        @test rd.tau == Float64.(tau_values(Octave, window, :adev))
        @test maximum(rd.tau) <= (window - 2) ÷ 2     # window clamp, ≪ record's
        @test all(isfinite, rd.dev)

        rh = dhdev(pd, Octave; window=window)
        @test rh.tau == Float64.(tau_values(Octave, window, :hdev))
        @test maximum(rh.tau) <= (window - 2) ÷ 3
        @test all(isfinite, rh.dev)
    end

    @testset "FrequencyData dispatch matches the cumsum-equivalent phase" begin
        Random.seed!(13)
        tau0 = 0.5
        y  = randn(1024) .* 1e-9
        fd = FrequencyData(y, tau0)
        pd = PhaseData(cumsum(y) .* tau0, tau0)
        ms = [1, 2, 4]

        rf = dadev(fd, ms; window=256)
        rp = dadev(pd, ms; window=256)
        @test rf.t == rp.t
        @test rf.tau == rp.tau == ms .* tau0
        @test rf.dev ≈ rp.dev
        @test dhdev(fd, ms; window=256).dev ≈ dhdev(pd, ms; window=256).dev
        # TauMode form on FrequencyData runs the same dispatch chain.
        @test dadev(fd, Octave; window=256).dev ≈ dadev(pd, Octave; window=256).dev
    end

    @testset "argument validation" begin
        pd = PhaseData(cumsum(randn(64)) .* 1e-9, 1.0)
        @test_throws ArgumentError dadev(pd, [1]; window=65)        # window > N
        @test_throws ArgumentError dadev(pd, [1]; window=3)         # below m=1 support
        @test_throws ArgumentError dhdev(pd, [1]; window=4)         # below m=1 support
        @test_throws ArgumentError dadev(pd, [1]; window=32, step=0)
        @test_throws ArgumentError dhdev(pd, [1]; window=32, step=-1)
    end

    @testset "compact show" begin
        pd = PhaseData(cumsum(randn(256)) .* 1e-9, 1.0)
        r = dadev(pd, [1, 2]; window=64)
        s = sprint(show, r)
        @test occursin("DynamicStabilityResult(:dadev", s)
        @test occursin("window=64", s)
        @test occursin("×2 map", s)
        # text/plain falls back to the same one-line summary (no table dump).
        @test sprint(show, MIME"text/plain"(), r) == s
    end
end
