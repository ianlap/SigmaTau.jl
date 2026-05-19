# models/ensemble.jl — Clock-Ensemble Time-Scale Model
#
# Stacked-state Kalman formulation of a clock ensemble per Stein 2003
# "Time Scales Demystified", §V — the joint state vector concatenates
# the per-clock states, Φ and Q are block-diagonal, and the (N−1)
# observable phase differences against a reference clock form H. The
# reference clock's absolute phase is non-observable (Stein §II); its
# covariance diagonal grows unbounded but every observable linear
# combination of states remains tight. This is the Galleani–Tavella
# free-clock formulation.
#
# The ensemble is itself an `AbstractClockModel`, so the existing
# `predict!` / `update!` loop on `StandardKalmanFilter` consumes it
# without modification — only `state_transition`, `process_noise`,
# `measurement_matrix`, and `measurement_noise` change shape.

"""
    EnsembleWeights{N}

Stein 2003 §VI–VII per-clock shock-allocation weights, normalised so
each weight vector sums to one. Stored on `ClockEnsemble` for use by
downstream shock-recombination steps (eqs. 6.2 / 6.3 / 7.2 / 7.3); the
joint stacked Kalman filter itself does not consume these — the
weights are an interpretation layer for recovering individual-clock
estimates from the joint state.

- `a::SVector{N,Float64}` — phase-shock (ε) weights, ∝ 1/q1.
- `b::SVector{N,Float64}` — frequency-shock (η) weights, ∝ 1/q2.
- `c::SVector{N,Float64}` — frequency-aging-shock (α) weights, ∝ 1/q3
  for clocks that model drift; identically zero for `TwoStateClock`
  ensembles and for `ThreeStateClock` ensembles whose members all have
  `q3 == 0`.
"""
struct EnsembleWeights{N}
    a::SVector{N,Float64}
    b::SVector{N,Float64}
    c::SVector{N,Float64}
end

# ── Stein auto-weights ───────────────────────────────────────────────────────

"""
    _stein_weights(clocks::NTuple{N, <:AbstractClockModel}) → EnsembleWeights{N}

Compute Stein §VI–VII inverse-noise weights from each clock's
diffusion coefficients (`a_i ∝ 1/q1_i`, `b_i ∝ 1/q2_i`,
`c_i ∝ 1/q3_i`) and normalise so each weight vector sums to one.
Throws `ArgumentError` when any clock has `q1 == 0` or `q2 == 0`
(Stein's inverse-variance weights are undefined in those cases — pass
explicit `weights` instead). `q3 == 0` is permitted: it collapses
the aging-weight vector to all zeros, the correct sentinel for a
no-drift ensemble.
"""
function _stein_weights(clocks::NTuple{N, ThreeStateClock}) where {N}
    any(c -> c.q1 == 0.0, clocks) &&
        throw(ArgumentError("ClockEnsemble auto-weights require q1 > 0 for every clock; pass explicit weights for the q1=0 case"))
    any(c -> c.q2 == 0.0, clocks) &&
        throw(ArgumentError("ClockEnsemble auto-weights require q2 > 0 for every clock; pass explicit weights for the q2=0 case"))

    a = SVector{N,Float64}(ntuple(i -> 1.0 / clocks[i].q1, N))
    b = SVector{N,Float64}(ntuple(i -> 1.0 / clocks[i].q2, N))
    c_raw = SVector{N,Float64}(ntuple(i -> clocks[i].q3 == 0.0 ? 0.0 : 1.0 / clocks[i].q3, N))
    c = sum(c_raw) == 0.0 ? zero(SVector{N,Float64}) : c_raw ./ sum(c_raw)

    return EnsembleWeights{N}(a ./ sum(a), b ./ sum(b), c)
end

function _stein_weights(clocks::NTuple{N, TwoStateClock}) where {N}
    any(c -> c.q1 == 0.0, clocks) &&
        throw(ArgumentError("ClockEnsemble auto-weights require q1 > 0 for every clock; pass explicit weights for the q1=0 case"))
    any(c -> c.q2 == 0.0, clocks) &&
        throw(ArgumentError("ClockEnsemble auto-weights require q2 > 0 for every clock; pass explicit weights for the q2=0 case"))

    a = SVector{N,Float64}(ntuple(i -> 1.0 / clocks[i].q1, N))
    b = SVector{N,Float64}(ntuple(i -> 1.0 / clocks[i].q2, N))

    return EnsembleWeights{N}(a ./ sum(a), b ./ sum(b), zero(SVector{N,Float64}))
end

# Fallback: any user-defined `AbstractClockModel` subtype hits this and
# gets a deliberate ArgumentError instead of a MethodError. Auto-weights
# are inverse-noise-coefficient specific to the polynomial-clock SDE;
# extend `_stein_weights` for new clock families before they can flow
# through `+` / auto-weights.
_stein_weights(::NTuple{N, <:AbstractClockModel}) where {N} =
    throw(ArgumentError("Stein auto-weights are only defined for TwoStateClock and ThreeStateClock ensembles. Pass explicit `weights::EnsembleWeights{N}` to ClockEnsemble, or extend `_stein_weights` for the new clock type."))

# ── ClockEnsemble type ───────────────────────────────────────────────────────

"""
    ClockEnsemble{N, M<:AbstractClockModel} <: AbstractClockModel

Homogeneous ensemble of `N` clocks for joint state-space Kalman
estimation under the Stein 2003 / Galleani–Tavella stacked
formulation. The joint state is the concatenation of each member's
state; Φ and Q are block-diagonal; H selects the `N-1` phase
differences against a reference clock.

Fields:
- `clocks::NTuple{N, M}` — the ensemble members (same concrete type).
- `weights::EnsembleWeights{N}` — Stein shock-allocation weights.
- `ref::Int` — reference-clock index (1-based; default 1). Phase
  differences are measured as `xᵢ − x_ref` for every `i ≠ ref`.
- `tau::Float64` — common discretisation step (every member must
  share this `tau`).

# Construction

```julia
ensemble = ClockEnsemble(clockA, clockB)                # 2-clock, ref=1, auto-weights
ensemble = ClockEnsemble((c1, c2, c3); ref=2)          # 3-clock, ref=2
ensemble = clockA + clockB                              # operator sugar
ensemble = clockA + clockB + clockC                     # left-associative
```

Heterogeneous mixes (e.g. `TwoStateClock + ThreeStateClock`) throw
`ArgumentError` — the joint state-space is only well-defined for
matching state dimensions. Mismatched `tau` between members likewise
throws.

# Use with `StandardKalmanFilter`

The ensemble is an `AbstractClockModel`, so the existing
`predict!` / `update!` / `prop!` API consumes it directly:

```julia
n   = nstates(ensemble)          # = N · nstates(M)
kf  = StandardKalmanFilter(zeros(n), Matrix(1e-12·I(n)))
for k in eachindex(z_diffs)
    predict!(kf, ensemble, ensemble.tau)
    update!(kf, ensemble, z_diffs[k])    # vector of length N-1 (or Real when N=2)
end
```

# References

S.R. Stein, "Time Scales Demystified", Proc. 2003 IFCS, pp. 223–227,
§V "Kalman-Filter Solution". L. Galleani and P. Tavella, "Time and
the Kalman filter", IEEE Control Systems Magazine, 2010 — joint
block-diagonal Φ/Q discussion.
"""
struct ClockEnsemble{N, M<:AbstractClockModel} <: AbstractClockModel
    clocks::NTuple{N, M}
    weights::EnsembleWeights{N}
    ref::Int
    tau::Float64

    function ClockEnsemble{N,M}(clocks::NTuple{N,M},
                                weights::EnsembleWeights{N},
                                ref::Int) where {N, M<:AbstractClockModel}
        N ≥ 2 || throw(ArgumentError("ClockEnsemble requires N ≥ 2 clocks; got N = $N"))
        1 ≤ ref ≤ N || throw(ArgumentError("ClockEnsemble ref index $ref out of range [1, $N]"))

        tau = clocks[1].tau
        for i in 2:N
            clocks[i].tau == tau ||
                throw(ArgumentError("ClockEnsemble members must share `tau`; got clocks[1].tau = $tau, clocks[$i].tau = $(clocks[i].tau)"))
        end
        return new{N,M}(clocks, weights, ref, tau)
    end
end

"""
    ClockEnsemble(clocks::NTuple{N,M}; ref=1, weights=nothing)

Primary constructor. When `weights === nothing`, auto-derives Stein
inverse-noise weights via [`_stein_weights`](@ref).
"""
function ClockEnsemble(clocks::NTuple{N,M};
                       ref::Int = 1,
                       weights::Union{Nothing,EnsembleWeights{N}} = nothing) where {N, M<:AbstractClockModel}
    w = weights === nothing ? _stein_weights(clocks) : weights
    return ClockEnsemble{N,M}(clocks, w, ref)
end

"""
    ClockEnsemble(c1::M, c2::M, more::M...; kwargs...)

Varargs sugar — builds the `NTuple` and delegates.
"""
function ClockEnsemble(c1::M, c2::M, more::M...; kwargs...) where {M<:AbstractClockModel}
    return ClockEnsemble((c1, c2, more...); kwargs...)
end

# ── Operator overloading: clockA + clockB → ClockEnsemble ────────────────────

"""
    Base.:+(c1::M, c2::M) where {M<:AbstractClockModel}

Build a two-clock `ClockEnsemble` (homogeneous, `ref=1`, auto-weights).
Chained sums extend the ensemble: `c1 + c2 + c3` produces a 3-clock
ensemble. Mixed types (`TwoStateClock + ThreeStateClock`) hit the
catch-all method below and throw.
"""
Base.:+(c1::M, c2::M) where {M<:AbstractClockModel} = ClockEnsemble((c1, c2))

Base.:+(e::ClockEnsemble{N,M}, c::M) where {N, M<:AbstractClockModel} =
    ClockEnsemble((e.clocks..., c); ref=e.ref)

Base.:+(c::M, e::ClockEnsemble{N,M}) where {N, M<:AbstractClockModel} =
    ClockEnsemble((c, e.clocks...); ref=e.ref + 1)

Base.:+(e1::ClockEnsemble{N,M}, e2::ClockEnsemble{K,M}) where {N, K, M<:AbstractClockModel} =
    ClockEnsemble((e1.clocks..., e2.clocks...); ref=e1.ref)

# Heterogeneous catch-all — fires when the typed methods above don't unify.
Base.:+(::AbstractClockModel, ::AbstractClockModel) =
    throw(ArgumentError("ClockEnsemble requires homogeneous clock types; mixing e.g. TwoStateClock and ThreeStateClock is not supported. Use a common model type or pre-convert."))

# ── State-space methods ──────────────────────────────────────────────────────

"""
    nstates(e::ClockEnsemble) → Int

Total joint-state dimension: `N · nstates(M)`.
"""
nstates(e::ClockEnsemble{N,M}) where {N,M} = N * nstates(e.clocks[1])

"""
    state_transition(e::ClockEnsemble, dt::Real) → SMatrix
    state_transition(e::ClockEnsemble)           → SMatrix

Block-diagonal stack of the per-clock transition matrices. Each
diagonal block is `state_transition(e.clocks[i], dt)`.
"""
function state_transition(e::ClockEnsemble{N,M}, dt::Real) where {N, M<:AbstractClockModel}
    ns   = nstates(e.clocks[1])
    Ntot = N * ns
    Φ    = MMatrix{Ntot, Ntot, Float64}(undef)
    fill!(Φ, 0.0)
    for k in 1:N
        block = state_transition(e.clocks[k], dt)
        off = (k - 1) * ns
        for i in 1:ns, j in 1:ns
            Φ[off + i, off + j] = block[i, j]
        end
    end
    return SMatrix(Φ)
end

state_transition(e::ClockEnsemble) = state_transition(e, e.tau)

"""
    process_noise(e::ClockEnsemble, dt::Real) → SMatrix
    process_noise(e::ClockEnsemble)           → SMatrix

Block-diagonal stack of the per-clock process-noise covariances.
Inter-clock cross-covariance is identically zero — the ensemble
assumes statistically independent clock noises, matching Stein §II
and the GTCH independence precondition used elsewhere in the package.
"""
function process_noise(e::ClockEnsemble{N,M}, dt::Real) where {N, M<:AbstractClockModel}
    ns   = nstates(e.clocks[1])
    Ntot = N * ns
    Q    = MMatrix{Ntot, Ntot, Float64}(undef)
    fill!(Q, 0.0)
    for k in 1:N
        block = process_noise(e.clocks[k], dt)
        off = (k - 1) * ns
        for i in 1:ns, j in 1:ns
            Q[off + i, off + j] = block[i, j]
        end
    end
    return SMatrix(Q)
end

process_noise(e::ClockEnsemble) = process_noise(e, e.tau)

"""
    measurement_matrix(e::ClockEnsemble) → SMatrix

Build the `(N-1) × (N·nstates(M))` phase-difference observation
matrix. For reference clock `ref`, each row corresponds to one
non-reference clock `i` and reads as
`zeros…  −1 at col 1+(ref-1)·ns  +1 at col 1+(i-1)·ns  zeros…` —
selecting `phase(i) − phase(ref)` from the joint state. Only the
phase component of each block is observed; the velocity/aging
components do not appear in the row, matching the phase-only
measurement convention of the underlying clock models.
"""
function measurement_matrix(e::ClockEnsemble{N,M}) where {N, M<:AbstractClockModel}
    ns   = nstates(e.clocks[1])
    Ntot = N * ns
    rows = N - 1
    H    = MMatrix{rows, Ntot, Float64}(undef)
    fill!(H, 0.0)

    ref_col = 1 + (e.ref - 1) * ns
    row = 0
    for i in 1:N
        i == e.ref && continue
        row += 1
        H[row, ref_col]               = -1.0
        H[row, 1 + (i - 1) * ns]      =  1.0
    end
    return SMatrix(H)
end

"""
    measurement_noise(e::ClockEnsemble) → SMatrix

`(N-1) × (N-1)` measurement-noise covariance for the phase
differences against the reference clock. Treating each clock's WPM
diffusion `q0_i` as the variance of an independent measurement
noise term `v_i`, every observation row is `z_i = (x_i − x_ref) +
(v_i − v_ref)`. The covariance of those independent noises gives:

- Diagonal: `R[i,i] = q0_ref + q0_i` — the independent-sum variance
  of one difference observation.
- Off-diagonal (`i ≠ j`): `R[i,j] = q0_ref` — every pair of
  differences shares the same `v_ref` term, so the cross-covariance
  is `Var(v_ref) = q0_ref`. Omitting this cross-term makes the filter
  overcount information from multiple differences against the same
  reference and underestimates uncertainty for `N ≥ 3`.

For `N = 2` the matrix is 1×1 and collapses to `q0_ref + q0_other`.
"""
function measurement_noise(e::ClockEnsemble{N,M}) where {N, M<:AbstractClockModel}
    rows = N - 1
    R    = MMatrix{rows, rows, Float64}(undef)
    q0_ref = e.clocks[e.ref].q0
    # Off-diagonals carry the shared reference-clock WPM term.
    fill!(R, q0_ref)
    # Diagonals add the non-reference clock's WPM.
    row = 0
    for i in 1:N
        i == e.ref && continue
        row += 1
        R[row, row] = q0_ref + e.clocks[i].q0
    end
    return SMatrix(R)
end
