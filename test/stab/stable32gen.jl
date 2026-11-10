# test/stab/stable32gen.jl — the reference kernels and bias-corrected TOTDEV
# against Stable32 and allantools on the 8192-point stable32gen.DAT record.
#
# Stable32 prints 4-5 significant figures, so the floor is rtol = 1e-4. It
# reports MTOTDEV uncorrected but bias-corrects TOTDEV (B = 1 − a·τ/T) and
# HTOTDEV (1/√(1+a)), so the raw kernels differ from it at the FLFM/RWFM rows
# by exactly those factors; `totdev` corrected for Stable32's noise type
# matches every row at 1e-4.
# allantools reports every estimator raw and is matched at 1e-11 (the CSV is
# written at %.17e), TOTDEV at 1e-7.

using Test
using SigmaTau

@testset "Stable32 on stable32gen.DAT" begin
    dir = joinpath(@__DIR__, "..", "fixtures", "validation")
    x = parse.(Float64, strip.(readlines(joinpath(dir, "stable32gen.DAT"))[11:end]))
    @test length(x) == 8192
    x_cs = pushfirst!(cumsum(x), 0.0)
    # Stable32 type => (raw deviation at m, rtol)
    raw = Dict(
        "Overlapping Allan"    => (m -> sqrt(LK.adev_var(x, m, 1.0)), 1e-4),
        "Modified Allan"       => (m -> sqrt(LK.mdev_var(x, m, 1.0, x_cs)), 1e-4),
        "Overlapping Hadamard" => (m -> sqrt(LK.hdev_var(x, m, 1.0)), 1e-4),
        "Time"                 => (m -> m * sqrt(LK.mdev_var(x, m, 1.0, x_cs)) / sqrt(3.0), 1e-4),
        "Modified Total"       => (m -> sqrt(LK.mtotdev_var(x, m, 1.0)), 5e-4),
        "Total"                => (m -> sqrt(LK.totdev_var_howe(x, m, 1.0)), 0.02),   # B ≤ 1.5 %
        "Hadamard Total"       => (m -> sqrt(LK.htotdev_var(x, m, 1.0)), 0.10),       # 1/√(1+a) ≤ 8.4 %
        "Allan"                => (m -> adev(PhaseData(x), [m]; ci=false, overlapping=false).dev[1], 1e-4),
        "Hadamard"             => (m -> hdev(PhaseData(x), [m]; ci=false, overlapping=false).dev[1], 1e-4))
    # CSV columns: Type, AF, Tau, N, Alpha, MinSigma, Sigma, MaxSigma
    rows = [(r[1], parse(Int, r[2]), parse(Float64, r[7]))
            for r in split.(readlines(joinpath(dir, "stable32out", "stable32_data_full.csv"))[2:end], ',')
            if length(r) >= 7 && haskey(raw, r[1])]
    @test length(rows) >= 50
    @test isempty([(k, m) for (k, m, ref) in rows
                   if !isapprox(raw[k][1](m), ref; rtol=raw[k][2])])

    pd = PhaseData(x, 1.0)
    # Non-overlapped ADEV / HDEV: Stable32's counts, and its intervals from the
    # stride-1 EDF at Stable32's noise type.
    csv = split.(readlines(joinpath(dir, "stable32out", "stable32_data_full.csv"))[2:end], ',')
    αsym = Dict(2 => :WHPM, 1 => :FLPM, 0 => :WHFM, -1 => :FLFM, -2 => :RWFM)
    for (k, d, est) in (("Allan", 2, adev), ("Hadamard", 3, hdev))
        s32 = filter(r -> r[1] == k, csv)
        ms = [parse(Int, r[2]) for r in s32]
        @test est(pd, ms; ci=false, overlapping=false).neff == [parse(Int, r[4]) for r in s32]
        edf = SigmaTau._nonoverlapped_edf(d, [αsym[parse(Int, r[5])] for r in s32], ms, length(x))
        lo, σ, hi = (parse.(Float64, getindex.(s32, c)) for c in 6:8)
        l, h = confidence_intervals(σ, edf, 0.683)
        @test all(isapprox.(l, lo; rtol=1e-3)) && all(isapprox.(h, hi; rtol=1e-3))
    end
    # From τ = 1024 Stable32 reads white FM, where σ(τ) runs flat and then
    # rises and SigmaTau reads flicker FM and random-walk FM.
    tot = filter(r -> r[1] == "Total", csv)
    ms = [parse(Int, r[2]) for r in tot]
    B = bias_correction([αsym[parse(Int, r[5])] for r in tot], :totvar, ms, length(x); published=true)
    @test length(tot) >= 10
    @test all(isapprox.(totdev(pd, ms; ci=false, correct_bias=false).dev ./ sqrt.(B),
                        parse.(Float64, getindex.(tot, 7)); rtol=1e-4))
    @test totdev(pd, [1024, 2048, 4096]).noise_type == [:FLFM, :RWFM, :RWFM]

    # correct_bias=:stable32 is Stable32's bias convention for each estimator:
    # its σ on the nine rows where the two programs read the same noise type.
    for (k, est) in (("Total", totdev), ("Modified Total", mtotdev), ("Time Total", ttotdev),
                     ("Hadamard Total", (p, m; kw...) -> htotdev(p, m; remove_drift=false, kw...)))
        s32 = filter(r -> r[1] == k, csv)
        r = est(pd, [parse(Int, r[2]) for r in s32]; correct_bias=:stable32)
        same = [αsym[parse(Int, r[5])] for r in s32] .== r.noise_type
        @test count(same) == 9
        @test r.dev[same] ≈ parse.(Float64, getindex.(s32, 7))[same] rtol = 1e-4
    end
    @test_throws ArgumentError totdev(pd, [4]; correct_bias=:sp1065)

    # allantools CSV columns: Type, AF, Tau, N, Sigma (raw totdev is _totdev_core)
    raw["Total"] = (m -> SigmaTau._totdev_core(x, [m], 1.0)[1], 0.0)
    at = [(r[1], parse(Int, r[2]), parse(Float64, r[5]))
          for r in split.(readlines(joinpath(dir, "allantools_out", "allantools_data_full.csv"))[2:end], ',')
          if length(r) >= 5 && r[5] != "nan" && haskey(raw, r[1])]
    @test length(at) >= 30
    @test isempty([(k, m) for (k, m, ref) in at
                   if !isapprox(raw[k][1](m), ref; rtol = k == "Total" ? 1e-7 : 1e-11)])
end
