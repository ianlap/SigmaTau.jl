using SigmaTau
using Random

# ── Load phase data ─────────────────────────────────────────────────────────
# From a Stable32-style .DAT file (10-line header, one value per line):
#
#   lines = readlines("my_data.DAT")
#   x     = parse.(Float64, strip.(lines[11:end]))
#   pd    = PhaseData(x, 1.0)          # 1.0 s sample interval
#
# Or from a plain vector:
Random.seed!(20260510)
x  = randn(512) .* 1e-9               # synthetic white-FM phase data
pd = PhaseData(x, 1.0)

# ── Compute Allan deviation ──────────────────────────────────────────────────
result = adev(pd, [1, 2, 4, 8, 16, 32])

println("τ (s):  ", result.tau)
println("ADEV:   ", result.dev)
println("noise:  ", result.noise_type)
println("CI low: ", result.ci_lower)
println("CI hi:  ", result.ci_upper)

# ── Overlay multiple deviations on one plot ──────────────────────────────────
# using Plots
# r_adev = adev(pd, [1, 2, 4, 8, 16, 32])
# r_mdev = mdev(pd, [1, 2, 4, 8, 16, 32])
# r_tdev = tdev(pd, [1, 2, 4, 8, 16, 32])
# plot(r_adev, label="ADEV")
# plot!(r_mdev, label="MDEV")
# plot!(r_tdev, label="TDEV")

# ── Save / load result ───────────────────────────────────────────────────────
save_result("my_adev.tsv", result)
result2 = load_result("my_adev.tsv")
println("\nRound-trip OK: ", result2.dev ≈ result.dev)
