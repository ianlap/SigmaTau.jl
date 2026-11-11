# Reference cross-check suite: every estimator, noise type α = +2 … −4, EDF,
# χ² confidence interval, noise ID, Welch spectrum, dynamic map, gap-tolerant
# ADEV and drift fits, against
#
#   * allantools (Python) on the committed records (fixtures/validation/
#     allantools/*.csv, regenerated in CI by tools/validation/gen_allantools.py),
#   * AllanDeviations.jl on the power-law records (fixtures/validation/
#     allandeviations/devs.csv, from tools/validation/gen_allandeviations.jl),
#   * Stable32 runs on PHASE.DAT (fixtures/validation/stable32_phasedat/) and
#     on a flicker PM frequency record (fixtures/validation/stable32_flickerpm/),
#   * the published values of the 1000-point test suite (NIST SP1065 §12.4),
#   * scipy.signal.welch for S_y / S_x and scipy.optimize.least_squares for
#     the drift fits.
#
# Each row carries the tolerance it is held to; rows with tol = Inf are
# documented convention differences that are reported (by
# tools/validation/report.jl) but not asserted. Point SIGMATAU_ALLANTOOLS_DIR
# (or SIGMATAU_ALLANDEVIATIONS_DIR) at a freshly generated directory to test
# against a live run.
using Test
using SigmaTau

isdefined(@__MODULE__, :CrossCheck) || include(joinpath(@__DIR__, "CrossCheck.jl"))
using .CrossCheck: Row, relerr

const _CC = CrossCheck

@testset "reference cross-check" begin
    for (label, f) in (
            ("allantools σ, every estimator × α = +2…−4", _CC.allantools_dev_rows),
            ("AllanDeviations.jl σ, every estimator × α", _CC.allandeviations_dev_rows),
            ("1000-point test suite, published values",  _CC.nbs14_rows),
            ("Stable32 PHASE.DAT σ, α and CI",            _CC.stable32_rows),
            ("Stable32 CI with Stable32's α",             _CC.stable32_ci_forced_alpha_rows),
            ("Stable32 σ on a flicker PM record",         _CC.stable32_flickerpm_rows),
            ("allantools Greenhall / SP1065 EDF",         _CC.allantools_edf_rows),
            ("allantools χ² interval factors",            _CC.allantools_ci_rows),
            ("allantools lag-1 ACF noise ID",             _CC.allantools_noiseid_rows),
            ("scipy Welch S_y / S_x",                     _CC.welch_rows),
            ("allantools per-window dynamic maps",        _CC.dynamic_rows),
            ("allantools gradev on records with gaps",    _CC.gap_rows),
            ("scipy least-squares drift fits",            _CC.drift_rows))
        @testset "$label" begin
            # One assertion per group; every row is still checked, and each
            # miss is logged with its reference value and tolerance.
            rows = filter(r -> !isinf(r.tol), f())
            misses = filter(r -> !_CC.passes(r), rows)
            for r in misses
                @info "cross-check miss" r.group r.dataset r.estimator r.m r.ref r.got relerr(r) r.tol
            end
            @test !isempty(rows)
            @test isempty(misses)
        end
    end
end
