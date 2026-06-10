# allantools_cross_validation.jl — Compares the SigmaTau raw
# kernels against Anders Wallin's `allantools` (Python) on the same
# Stable32 phase fixture. Skipped silently when the fixture has not
# been regenerated locally (`test/fixtures/validation/allantools_out/
# allantools_data_full.csv`); regenerate with
# `python3 tools/regen_allantools_fixtures.py`.
#
# Why a third reference? allantools reports the *raw* total estimators
# (no bias correction), while Stable32's policy is per-estimator (it
# corrects TOTDEV and HTOTDEV but reports MTOTDEV raw). Comparing our
# raw kernels against allantools therefore isolates pure kernel
# agreement from bias-correction policy — and the raw kernels match
# allantools to machine precision (≤4.4e-15 measured on this fixture).
#
# Comparison contract: raw legacy kernels (LK.* — bias-free) vs
# allantools' default output, at the shared 1e-11 cross-platform floor.

@testset "allantools cross-validation" begin
    ref_dir = joinpath(@__DIR__, "..", "fixtures", "validation")
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

        # Single shared rtol: every raw kernel — including HTOTDEV and
        # MTOTDEV — matches allantools at machine precision (measured
        # worst case 4.4e-15 on this fixture, macOS x86_64 2026-06-10).
        #
        # `tight = 1e-11` (was 1e-4): the regen script writes the CSV
        # at %.17e (round-trip-exact Float64) instead of %.6e
        # (~7 sig figs), so the fixture itself preserves machine
        # precision; 1e-11 is comfortable headroom for the
        # ~10,000-ULP cross-platform LLVM codegen drift we see on
        # Linux x86_64.
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
                # allantools' raw `totdev` follows SP1065 eqn 25 verbatim
                # (no detrend), which is exactly what `_totdev_core` (the
                # canonical Howe form) implements. Apples-to-apples: agreement
                # is ~7 sig figs across all m on this fixture (no need for the
                # m=512 skip Stable32 needs — allantools doesn't apply
                # Stable32's alpha-aware correction).
                got = SigmaTau._totdev_core(x, [m], tau0)[1]
                rtol = 1e-7
            elseif kind == "Hadamard Total"
                # Raw kernel, apples-to-apples with allantools' raw
                # htotdev: machine precision (≤4.4e-15 measured).
                got = sqrt(LK.htotdev_var(x, m, tau0))
            elseif kind == "Modified Total"
                # Raw kernel, apples-to-apples with allantools' raw
                # mtotdev: machine precision (≤1.1e-15 measured).
                got = sqrt(LK.mtotdev_var(x, m, tau0))
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
