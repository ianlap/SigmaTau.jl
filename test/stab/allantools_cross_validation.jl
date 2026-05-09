# allantools_cross_validation.jl — Compares the SigmaTau.Stab raw
# kernels against Anders Wallin's `allantools` (Python) on the same
# Stable32 phase fixture. Skipped silently when the fixture has not
# been regenerated locally (`reference/validation/allantools_out/
# allantools_data_full.csv`); regenerate with
# `python3 tools/regen_allantools_fixtures.py`.
#
# Why a third reference? Stable32 reports unbiased totals while our API
# applies the SP1065 bias factor; allantools (default settings) sits
# closer to Stable32 for raw-kernel comparison and lets us isolate
# bias-policy disagreement from boundary-policy disagreement.
#
# Comparison contract: raw legacy kernels (LK.* — bias-free) vs
# allantools' default output. Tolerances mirror the existing Stable32
# testset (`runtests.jl`), since the policy disagreements between
# allantools and Stable32 on totals are documented and small.

@testset "allantools cross-validation" begin
    ref_dir = joinpath(@__DIR__, "..", "..", "..", "reference", "validation")
    dat_path = joinpath(ref_dir, "stable32gen.DAT")
    at_csv   = joinpath(ref_dir, "allantools_out", "allantools_data_full.csv")

    if !isfile(at_csv)
        @info "allantools fixture not present (run tools/regen_allantools_fixtures.py); skipping"
    elseif !isfile(dat_path)
        @info "Stable32 phase fixture not present; skipping allantools cross-check"
    else
        # Phase fixture (10-line header, 8192 samples) — same as the
        # Stable32 testset.
        lines = readlines(dat_path)
        x = parse.(Float64, strip.(lines[11:end]))
        @test length(x) == 8192
        tau0 = 1.0
        x_cs = pushfirst!(cumsum(x), 0.0)

        # Allantools CSV columns: Type, AF, Tau, N, Sigma.
        rows = [split(line, ',') for line in readlines(at_csv)[2:end]]

        # Per-kernel rtol. ADEV/MDEV/HDEV/TDEV agree tightly between
        # SigmaTau and allantools (same kernel definition, same
        # boundary handling). TOTDEV / HTOTDEV / MTOTDEV agree less
        # tightly — different boundary-extension conventions.
        #
        # `tight = 1e-11` (was 1e-4): the regen script now writes the
        # CSV at %.17e (round-trip-exact Float64) instead of %.6e
        # (~7 sig figs), so the fixture itself preserves machine
        # precision. Three-way verification on macOS x86_64 (2026-05-08)
        # shows ours/legacy/allantools agree to ≤ 8.5e-14 worst case
        # on this fixture, so 1e-11 is comfortable headroom for the
        # ~10,000-ULP cross-platform LLVM codegen drift we see on
        # Linux x86_64. TOTDEV/HTOTDEV/MTOTDEV stay at their original
        # boundary-policy floors below.
        tight = 1e-11

        n_checked = Dict{String,Int}()
        n_skipped = 0
        for row in rows
            length(row) < 5 && continue
            kind = String(row[1])
            m    = parse(Int, row[2])
            sigma_ref_str = String(row[5])

            if sigma_ref_str == "nan"
                n_skipped += 1
                continue
            end
            sigma_ref = parse(Float64, sigma_ref_str)

            got = NaN
            rtol = tight
            if kind == "Overlapping Allan"
                got = sqrt(LK.adev_var(x, m, tau0))
            elseif kind == "Modified Allan"
                got = sqrt(LK.mdev_var(x, m, tau0, x_cs))
            elseif kind == "Overlapping Hadamard"
                got = sqrt(LK.hdev_var(x, m, tau0))
            elseif kind == "Time"
                # TDEV = τ · MDEV / √3
                mdev_v = sqrt(LK.mdev_var(x, m, tau0, x_cs))
                got = (m * tau0) * mdev_v / sqrt(3.0)
            elseif kind == "Total"
                got = sqrt(LK.totdev_var(x, m, tau0))
                rtol = 0.15   # boundary-extension policy floor (matches Stable32 testset)
            elseif kind == "Hadamard Total"
                got = sqrt(LK.htotdev_var(x, m, tau0))
                rtol = 0.10   # ~0.5% bias + boundary effects
            elseif kind == "Modified Total"
                got = sqrt(LK.mtotdev_var(x, m, tau0))
                rtol = 0.05   # raw-kernel match per comparison_report.md
            else
                continue
            end

            @test got ≈ sigma_ref rtol=rtol
            n_checked[kind] = get(n_checked, kind, 0) + 1
        end

        @test sum(values(n_checked)) >= 30
        @info "allantools cross-validation: " *
              join(["$k=$v" for (k, v) in sort(collect(n_checked); by=first)], ", ") *
              (n_skipped > 0 ? " (skipped $n_skipped NaN rows)" : "")
    end
end
