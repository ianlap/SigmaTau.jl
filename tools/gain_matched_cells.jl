# tools/gain_matched_cells.jl — matched-cell EDF gain comparison backing the
# paper's gain claims: MHTOT/MHDEV vs MTOT/MDEV, with every ratio formed from
# measured cells on the SAME realization stream — no fitted EDF model enters
# anywhere. Fitted-model gain ratios were shown (2026-08-13) to manufacture
# spurious sub-unity mid-range readings on both sides of the comparison.
#
# Usage: julia --project=. --threads=auto tools/gain_matched_cells.jl
#
# α ∈ {2, 0, -1, -2}, N = 4096, R = 2000, octave m plus edge cells. Bump R for
# tighter cells (SE ≈ √(2·(2+12/ν)/R) per gain entry).
using Random, Statistics, FFTW, SigmaTau
using SigmaTau: _gen_powerlaw_y

const N = 4096
const R = 2000
const ms = [32, 64, 128, 256, 512, 640, 768, 896]

# Process-stable seed derivation (splitmix64 finalizer), as in mc_mhtotdev.jl.
function _seed(parts::Integer...)
    h = 0x243f6a8885a308d3
    for p in parts
        h ⊻= reinterpret(UInt64, Int64(p))
        h = (h ⊻ (h >> 30)) * 0xbf58476d1ce4e5b9
        h = (h ⊻ (h >> 27)) * 0x94d049bb133111eb
        h ⊻= (h >> 31)
    end
    return h
end

nu(v) = 2 * mean(v)^2 / var(v)

for alpha in (2, 0, -1, -2)
    nc = length(ms)
    VH = Matrix{Float64}(undef, nc, R); WH = Matrix{Float64}(undef, nc, R)
    VM = Matrix{Float64}(undef, nc, R); WM = Matrix{Float64}(undef, nc, R)
    Threads.@threads :dynamic for r in 1:R
        rng = Xoshiro(_seed(8801, alpha, N, r))
        pd = PhaseData(cumsum(_gen_powerlaw_y(float(alpha), N; rng=rng)), 1.0)
        VH[:, r] = mhtotdev(pd, ms; ci=false, correct_bias=false).dev .^ 2
        WH[:, r] = mhdev(pd, ms; ci=false).dev .^ 2
        VM[:, r] = mtotdev(pd, ms; ci=false, correct_bias=false).dev .^ 2
        WM[:, r] = mdev(pd, ms; ci=false).dev .^ 2
    end
    println("alpha=$alpha  N=$N  R=$R")
    println("m     T/tau  gain_MHTOT/MHDEV  gain_MTOT/MDEV   se%")
    for k in eachindex(ms)
        m = ms[k]; tt = (N - 1) / m
        gh = nu(@view VH[k, :]) / nu(@view WH[k, :])
        gm = nu(@view VM[k, :]) / nu(@view WM[k, :])
        nuh = nu(@view VH[k, :])
        se = 100 * sqrt(2 * (2 + 12 / nuh) / R)
        println(rpad(m, 6), rpad(round(tt, digits=1), 7),
                rpad(round(gh, digits=3), 18), rpad(round(gm, digits=3), 16),
                round(se, digits=1))
    end
end
