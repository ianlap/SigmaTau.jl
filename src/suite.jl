# `stability`: a suite of deviations in one call, and `ncornered_hat`, the N-cornered hat.
"""
Default deviation set computed by [`stability`](@ref) when `devs` is not given:
overlapping Allan, modified Allan, overlapping Hadamard, and modified Hadamard.
All four are fractional-frequency (σ_y) quantities, so the default suite
overlays on one ordinate; time deviations (`tdev`, `htdev`) are σ_x quantities
in seconds and must be requested explicitly.
"""
const DEFAULT_DEVIATIONS = (:adev, :mdev, :hdev, :mhdev)

function _stability(data, N::Int, kind::Symbol, tau0::Float64;
                    devs = DEFAULT_DEVIATIONS, m = Octave, taus = nothing,
                    ci::Bool = true, confidence::Real = DEFAULT_CONFIDENCE,
                    correct_bias::Bool = true, remove_drift::Bool = true)
    if taus !== nothing
        Base.depwarn("stability: `taus` is deprecated, use `m`", :stability)
        m = taus
    end
    results = StabilityResult[]
    for sym in devs
        name = Symbol(sym)
        ms = m isa TauMode ? _resolve_m(name, m, N) : collect(Int, m)
        push!(results, deviation(name, data, ms; ci, confidence, correct_bias, remove_drift))
    end
    tau_mode = m isa TauMode ? Symbol(string(m)) : :explicit
    return StabilitySuite(results, kind, tau0, N, ci ? Float64(confidence) : NaN, tau_mode)
end

"""
    stability(data::PhaseData; devs=DEFAULT_DEVIATIONS, m=Octave,
              ci=true, confidence=DEFAULT_CONFIDENCE, correct_bias=true,
              remove_drift=true) → StabilitySuite

Compute a suite of stability deviations from one record in a single call, the
batch analog of `adev`, `mdev`, …. Each symbol in `devs` is run through
[`deviation`](@ref) with a per-deviation averaging-factor grid derived from
`m`, so every deviation gets a grid valid for its own algorithmic m-max.
Returns a [`StabilitySuite`](@ref). Also accepts `FrequencyData`.

`devs` is any iterable of estimator names from [`estimators`](@ref)
(symbols or strings). `m` is a [`TauMode`](@ref) (e.g. `Octave`, `Decade`)
or an integer vector of averaging factors applied to every deviation, as in
the estimator functions.
`correct_bias` reaches only the estimators that take it (the total family and
`theo1`), and `remove_drift` only `mhtotdev`. Note that an explicit `m`
vector is passed to every deviation verbatim: `:theo1` interprets it as even
Thêo1 factors (odd entries give NaN rows) and `:theoh` as τ-grid factors, so
prefer a `TauMode` when mixing the Theo family with others.

# Examples

```jldoctest
julia> using SigmaTau

julia> pd = PhaseData(collect(1.0:64.0), 1.0);

julia> suite = stability(pd; devs=(:adev, :hdev), m=Octave, ci=false);

julia> keys(suite)
2-element Vector{Symbol}:
 :adev
 :hdev

julia> suite[:adev].deviation_type
:adev
```
"""
stability(data::PhaseData; kwargs...) =
    _stability(data, length(data.x), :phase, data.tau0; kwargs...)

stability(data::FrequencyData; kwargs...) =
    _stability(data, length(data.y) + 1, :frequency, data.tau0; kwargs...)   # N + 1 phase points


# N-cornered hat

# Resolve the measured pair (i, j), i < j, from a possibly part-filled
# matrix: prefer the upper-triangle entry, accept the transposed one
# (s²(x_i − x_j) = s²(x_j − x_i), so orientation is irrelevant).
function _nch_pair(pairwise::AbstractMatrix{StabilityResult}, i::Int, j::Int)
    isassigned(pairwise, i, j) && return pairwise[i, j]
    isassigned(pairwise, j, i) && return pairwise[j, i]
    throw(ArgumentError("ncornered_hat: no pairwise result for clocks ($i, $j); fill pairwise[$i, $j] (or its transpose)"))
end

"""
    ncornered_hat(pairwise::AbstractMatrix{StabilityResult}) → Vector{StabilityResult}
    ncornered_hat(r_ab, r_bc, r_ca) → Vector{StabilityResult}

N-cornered-hat noise separation (the classic Gray–Allan three-cornered
hat at `N = 3`; see NIST SP1065 §10.7 [riley-2008-sp1065](@cite)).
Given the pairwise deviations between every pair of `N ≥ 3` clocks,
estimates each clock's *individual* deviation under the assumption that
the clock noises are statistically independent, so the measured pairwise
variances satisfy `s²ᵢⱼ(τ) = σ²ᵢ(τ) + σ²ⱼ(τ)`.

With `Rᵢ = Σ_{j≠i} s²ᵢⱼ` (clock `i`'s row sum) and `S = Σ_{j<k} s²ⱼₖ`
(the sum over all pairs), the separation is

```math
\\hat\\sigma^2_i(\\tau) = \\frac{1}{N-2}\\left( R_i - \\frac{S}{N-1} \\right)
```

At `N = 3` this is the classic three-cornered hat,
`σ̂²_A = ½(s²_AB + s²_AC − s²_BC)`.

`pairwise[i, j]` (for `i < j`; the transposed entry is also accepted, so
fill either triangle, or both) is the deviation computed on the
difference record between clocks `i` and `j`. All entries must share the
same `deviation_type` and the same τ grid; otherwise an `ArgumentError`
is thrown. The 3-argument convenience form takes the three pairwise
results in the cyclic order A−B, B−C, C−A and returns `[A, B, C]`
(the orientation of each difference record does not matter).

Returned conventions (one `StabilityResult` per clock, in clock order):

- `deviation_type` is kept from the inputs: the statistic is still an
  ADEV, say; the hat only separates the pairwise variances into per-clock
  contributions.
- Non-positive variance estimates (common when one clock dominates or
  the independence assumption fails) are reported as `NaN` in `dev`
  rather than clamped to zero; they are invalid estimates, and zeros
  break log-scale plots.
- `ci`, `edf`, and `noise_type` are empty: there is no established
  confidence model for the hat-separated estimates (the inversion mixes
  χ² variables of differing sign), so no bounds are fabricated.
- `neff` is the elementwise minimum over all pairwise inputs: every
  separated estimate uses every measured pair, so the least-averaged
  input bounds its quality.

# Examples

```jldoctest
julia> using SigmaTau

julia> x1, x2, x3 = randn(256), randn(256), randn(256);

julia> r(a, b) = adev(PhaseData(a .- b, 1.0), [1, 2, 4]; ci=false);

julia> clocks = ncornered_hat(r(x1, x2), r(x2, x3), r(x3, x1));

julia> length(clocks), clocks[1].deviation_type
(3, :adev)
```
"""
function ncornered_hat(pairwise::AbstractMatrix{StabilityResult})
    n_clocks, n_cols = size(pairwise)
    n_clocks == n_cols || throw(ArgumentError(
        "ncornered_hat: pairwise matrix must be square, got $(n_clocks)×$(n_cols)"))
    n_clocks >= 3 || throw(ArgumentError(
        "ncornered_hat: need at least 3 clocks for a cornered hat, got $n_clocks"))

    ref = _nch_pair(pairwise, 1, 2)
    n_tau = length(ref.tau)

    # Accumulate per-τ row sums Rᵢ and the all-pairs sum S, validating
    # every pair against the reference grid and type.
    R = zeros(n_clocks, n_tau)
    S = zeros(n_tau)
    neff_min = fill(typemax(Int), n_tau)
    for i in 1:(n_clocks - 1), j in (i + 1):n_clocks
        r = _nch_pair(pairwise, i, j)
        r.deviation_type === ref.deviation_type || throw(ArgumentError(
            "ncornered_hat: pairwise ($i, $j) is :$(r.deviation_type) but (1, 2) is :$(ref.deviation_type); all pairs must share one deviation type"))
        r.tau == ref.tau || throw(ArgumentError(
            "ncornered_hat: pairwise ($i, $j) was computed on a different τ grid; all pairs must share one grid"))
        @inbounds for k in 1:n_tau
            v = r.dev[k]^2
            R[i, k] += v
            R[j, k] += v
            S[k] += v
            neff_min[k] = min(neff_min[k], r.neff[k])
        end
    end

    # σ̂²ᵢ = (Rᵢ − S/(N−1)) / (N−2); non-positive estimates become NaN.
    results = Vector{StabilityResult}(undef, n_clocks)
    for i in 1:n_clocks
        dev_i = Vector{Float64}(undef, n_tau)
        for k in 1:n_tau
            v = (R[i, k] - S[k] / (n_clocks - 1)) / (n_clocks - 2)
            dev_i[k] = v > 0 ? sqrt(v) : NaN
        end
        results[i] = StabilityResult(ref.deviation_type, copy(ref.tau), dev_i,
                                     Symbol[], _CIBound[], Float64[], copy(neff_min),
                                     ref.tau0, ref.n, NaN)
    end
    return results
end

# Three-cornered hat from A−B, B−C, C−A; C−A fills pairwise[1, 3] since
# orientation is irrelevant.
function ncornered_hat(r_ab::StabilityResult, r_bc::StabilityResult, r_ca::StabilityResult)
    pairwise = Matrix{StabilityResult}(undef, 3, 3)
    pairwise[1, 2] = r_ab
    pairwise[2, 3] = r_bc
    pairwise[1, 3] = r_ca
    return ncornered_hat(pairwise)
end
