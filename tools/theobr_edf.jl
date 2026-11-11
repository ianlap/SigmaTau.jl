# tools/theobr_edf.jl — scatter of the ThêoBR ratio on simulated FM noise and
# its correlation with Thêo1: the source of `_THEOBR_SCATTER` in src/edf.jl,
# and a check of `_theobr_edf` built on it.
#
#   julia -t auto --project=. tools/theobr_edf.jl        # about two minutes
#
# ThêoBR is R·Thêo1(m), with R the mean ratio of `_theobr_ratio` over the
# ladder to T/2 or over Howe's ladder to T/10. Part 1 prints the relative
# scatter c of R by record length, beside c∞·√(1 + (N0/N)^1.5). Part 2 prints
# ρ, the correlation of ln R with ln Thêo1(m), at τ/T = `_THEOBR_U`. The
# table in src/edf.jl is part 2 rounded, with zero at τ/T = 0 and for the
# small negative values; c∞ and N0 were set by eye on part 1.
# Part 3 compares the model with the EDF measured on other records,
# 2·mean²/var of the corrected variance, and counts how often a 68.3 %
# interval built on the model holds the true Allan variance.
using SigmaTau, Distributions, Printf, Random, Statistics
using SigmaTau: _theobr_ratio, _theo1_core, _theobr_edf, _THEOBR_U, _THEOBR_SCATTER

const NAMES = ("white FM", "flicker FM", "random-walk FM")
const LADDERS = (("T/2", :half, false), ("T/10", :tenth, true))

# R on both ladders and the Thêo1 variance at `ms`, for `nrec` records.
function sample(α, N, ms, nrec, seed)
    R = zeros(nrec, 2); T = zeros(nrec, length(ms))
    Threads.@threads for s in 1:nrec
        x = Vector{Float64}(noise_gen(PhaseData, N, 1.0; h=Dict(α => 1.0), rng=Xoshiro(seed + s)).x)
        R[s, 1] = _theobr_ratio(x, 1.0); R[s, 2] = _theobr_ratio(x, 1.0; howe=true)
        T[s, :] = _theo1_core(x, ms, 1.0) .^ 2
    end
    return R, T
end
factors(N, us) = [clamp(2 * round(Int, u * N / 1.5), 10, 2 * ((N - 1) ÷ 2)) for u in us]

println("1. Relative scatter of R, measured and modeled")
for (k, α) in enumerate((0, -1, -2)), (j, (name, key, _)) in enumerate(LADDERS)
    s = _THEOBR_SCATTER[key]
    @printf("  %-14s %-4s", NAMES[k], name)
    for (N, nrec) in ((64, 4000), (128, 4000), (256, 4000), (512, 4000), (2048, 2000), (8192, 400))
        R, _ = sample(α, N, [10], nrec, 1000N + 100k)
        @printf("  N = %d: %.3f %.3f", N, std(R[:, j]) / mean(R[:, j]), s.c[k] * sqrt(1 + (s.n0 / N)^1.5))
    end
    println()
end

println("2. Correlation of ln R with ln Thêo1 at τ/T = ", _THEOBR_U[2:end], ", N = 2048")
for (k, α) in enumerate((0, -1, -2))
    R, T = sample(α, 2048, factors(2048, _THEOBR_U[2:end]), 4000, 5000k)
    for (j, (name, _, _)) in enumerate(LADDERS)
        @printf("  %-14s %-4s", NAMES[k], name)
        for i in axes(T, 2); @printf(" %5.2f", cor(log.(R[:, j]), log.(T[:, i]))); end
        println()
    end
end

println("3. EDF measured and modeled, and coverage of the 68.3 % interval, N = 300 and 1024")
covers(v, truth, ν) = v * ν / quantile(Chisq(ν), 0.8415) <= truth <= v * ν / quantile(Chisq(ν), 0.1585)
for N in (300, 1024), (k, α) in enumerate((0, -1, -2))
    us = (0.02, 0.1, 0.2, 0.4, 0.6); ms = factors(N, us)
    R, T = sample(α, N, ms, 4000, 9000N + 100k)
    # True Allan variance at τ = 0.75m: 1/(2τ) and (4π²/6)·τ·(1 + 1/(2τ²)) for
    # white and random-walk FM with h = 1; a level for flicker FM, taken from
    # AVAR at τ = N/32 … N/8 on the same records' generator.
    flfm = α == -1 ? mean(mean(adev(noise_gen(PhaseData, N, 1.0; h=Dict(-1 => 1.0), rng=Xoshiro(77N + s)), [N ÷ 32, N ÷ 16, N ÷ 8]; ci=false).dev .^ 2) for s in 1:2000) : 0.0
    truth(τ) = α == 0 ? 1 / 2τ : α == -1 ? flfm : (2π^2 / 3) * τ * (1 + 1 / (2τ^2))
    for (j, (name, _, howe)) in enumerate(LADDERS)
        @printf("  N = %-4d %-14s %-4s", N, NAMES[k], name)
        for (i, m) in enumerate(ms)
            V = R[:, j] .* T[:, i]; ν = _theobr_edf(α, m, N; howe=howe)
            @printf("  %.2f T: %5.1f %5.1f, %2.0f %%", us[i], 2mean(V)^2 / var(V), ν, 100mean(covers(v, truth(0.75m), ν) for v in V))
        end
        println()
    end
end
