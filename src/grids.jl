# grids.jl — averaging-factor grids and conversion helpers.
#
# Merged from stab/taus.jl (TauMode, _kernel_m_max, _grid, tau_values) and
# stab/utils.jl (_f64, _freq_to_phase, _phase_to_freq, _default_m_values).


# ──────────────────────────────────────────────────────────────────────
# ── stab/taus.jl ────────────────────────────────────────────────────────

# stab/taus.jl — Averaging-factor grid selection (`TauMode` + `tau_values`)

"""
    TauMode

Selector for the averaging-factor grid used by the public deviation API and
[`stability`](@ref). Pass an instance in place of an explicit `m_values` vector,
e.g. `adev(pd, Octave)`.

Instances:

| `TauMode`       | Spacing            | Example grid (large `m_max`)        |
|-----------------|--------------------|-------------------------------------|
| `AllTaus`       | every integer      | `1, 2, 3, 4, 5, …`                  |
| `Octave`        | factor 2           | `1, 2, 4, 8, 16, …`                 |
| `HalfOctave`    | factor √2          | `1, 2, 3, 4, 6, 8, 11, 16, …`       |
| `QuarterOctave` | factor 2^(1/4)     | `1, 2, 3, 4, 5, 6, 7, 8, 10, …`     |
| `Decade`        | factor 10          | `1, 10, 100, 1000, …`               |
| `HalfDecade`    | factor √10         | `1, 3, 10, 32, 100, …`              |

`Octave` is the package default and reproduces the grid that the zero-arg
convenience methods (`adev(pd)`, …) have always used. Grids are rounded to
integers, deduplicated, and clamped to each kernel's algorithmic m-max (see
[`tau_values`](@ref)).
"""
@enum TauMode AllTaus Octave HalfOctave QuarterOctave Decade HalfDecade

"""
    _kernel_m_max(N::Int, kernel::Symbol) → Int

Largest averaging factor `m` for which `kernel`'s core L-check still yields at
least one analysis window, given `N` samples. Shared by [`tau_values`](@ref)
and [`_default_m_values`](@ref).

Per-kernel m-max (derived from the `L`/`Ne` guard in each `_*_core`):

| `kernel`                  | m_max          |
|---------------------------|----------------|
| `:adev`, `:pdev`          | `(N − 2) ÷ 2`  |
| `:totdev`                 | `(N − 1) ÷ 2`  |
| `:mdev`, `:tdev`          | `(N − 1) ÷ 3`  |
| `:mtotdev`, `:ttotdev`    | `N ÷ 3`        |
| `:hdev`                   | `(N − 2) ÷ 3`  |
| `:htotdev`                | `(N − 1) ÷ 3`  |
| `:mhdev`, `:htdev`        | `(N − 1) ÷ 4`  |
| `:mhtotdev`               | `N ÷ 4`        |
| `:mtie`                   | `N − 1`        |

The ordinary and modified estimators require ≥2 analysis windows at the largest
τ — a single window is one raw difference (EDF ≈ 1) — so their caps sit one step
below the total family, whose subsequence extension keeps a single subsequence a
valid long-τ estimate. (HTOTDEV's general branch runs on `y = diff(x)` of length
`N−1`, giving `(N−1) ÷ 3`; MTOTDEV runs on phase directly, giving `N ÷ 3`.)

Throws `ArgumentError` for an unknown kernel symbol or an `N` too short to admit
any `m ≥ 1`.
"""
function _kernel_m_max(N::Int, kernel::Symbol)
    # Ordinary/modified estimators use raw differences, so the largest τ must
    # retain ≥2 analysis windows — a single window is one difference (EDF ≈ 1)
    # and is rejected by the cores. Total estimators extend each subsequence
    # (reflection / time reversal), so a single subsequence is still a valid
    # long-τ estimate (the point of the total approach); they keep the natural
    # ≥1-subsequence reach.
    m_max = if kernel === :adev || kernel === :pdev
        (N - 2) ÷ 2
    elseif kernel === :totdev
        (N - 1) ÷ 2
    elseif kernel === :mdev || kernel === :tdev
        (N - 1) ÷ 3
    elseif kernel === :mtotdev || kernel === :ttotdev
        N ÷ 3
    elseif kernel === :hdev
        (N - 2) ÷ 3
    elseif kernel === :htotdev
        (N - 1) ÷ 3
    elseif kernel === :mhdev || kernel === :htdev
        (N - 1) ÷ 4
    elseif kernel === :mhtotdev
        N ÷ 4
    elseif kernel === :mtie
        N - 1
    else
        throw(ArgumentError("tau_values: unknown kernel symbol :$kernel"))
    end
    m_max < 1 && throw(ArgumentError(
        "tau_values: N=$N is too short to support any m for :$kernel"))
    return m_max
end

"""
    _grid(mode::TauMode, m_max::Int) → Vector{Int}

Averaging-factor grid for `mode`, clamped to `[1, m_max]`. `Octave` is returned
via the exact integer formula `[2^k for k in 0:floor(Int, log2(m_max))]` so the
default grid is byte-identical to the historical one; the other geometric modes
round their factor powers to `Int`, deduplicate, and clamp.
"""
function _grid(mode::TauMode, m_max::Int)
    if mode === AllTaus
        return collect(1:m_max)
    elseif mode === Octave
        return [2^k for k in 0:floor(Int, log2(m_max))]
    end
    f = mode === HalfOctave    ? sqrt(2.0) :
        mode === QuarterOctave ? 2.0^(1 / 4) :
        mode === Decade        ? 10.0 :
        sqrt(10.0)   # HalfDecade
    ms = Int[]
    k = 0
    while (m = max(1, round(Int, f^k))) <= m_max
        push!(ms, m)
        k += 1
    end
    return unique(ms)
end

"""
    tau_values(mode::TauMode, N::Int, kernel::Symbol) → Vector{Int}

Averaging factors `m` (with `τ = m·τ₀`) for the grid `mode`, bounded above by
`kernel`'s algorithmic m-max for an `N`-sample record. The return value is the
same `Vector{Int}` accepted by the explicit `m_values` form of every deviation
(e.g. `adev(pd, tau_values(Decade, length(pd.x), :adev))`), and the deviation
functions accept a `TauMode` directly (`adev(pd, Decade)`) as shorthand.

Note the name reads like seconds, but the returned values are averaging factors,
matching the existing `m_values` convention. `kernel` is the deviation symbol
(`:adev`, `:mdev`, …) that determines the m-max clamp; see [`_kernel_m_max`](@ref).

`tau_values(Octave, N, kernel)` equals [`_default_m_values`](@ref)`(N, kernel)`.
"""
tau_values(mode::TauMode, N::Int, kernel::Symbol) = _grid(mode, _kernel_m_max(N, kernel))


# ──────────────────────────────────────────────────────────────────────
# ── stab/utils.jl ───────────────────────────────────────────────────────

# utils.jl — Shared helpers for the stability API

"""
    _f64(v::AbstractVector{<:Real}) → Vector{Float64}

Promote a sample vector to the `Vector{Float64}` the core kernels and
`identify_noise` require. A no-op (returns the same array, no copy) when `v`
is already a `Vector{Float64}`; otherwise allocates a converted copy. Lets
`PhaseData{Float32}` (and any other `AbstractFloat` element type) flow through
the public deviation API instead of hitting a `MethodError` on the
`Vector{Float64}`-typed kernels.
"""
_f64(v::Vector{Float64}) = v
_f64(v::AbstractVector{<:Real}) = Vector{Float64}(v)

"""
    _freq_to_phase(data::FrequencyData) → PhaseData

Convert fractional-frequency samples to phase residuals using the running
integral `x[k] = τ₀ · Σⱼ₌₁ᵏ y[j]`. Length is preserved (N → N), matching the
convention that `adev(FrequencyData(y, τ₀), …)` agrees with
`adev(PhaseData(cumsum(y)·τ₀, τ₀), …)`.
"""
_freq_to_phase(data::FrequencyData) = PhaseData(cumsum(data.y) .* data.tau0, data.tau0)

"""
    _phase_to_freq(data::PhaseData) → FrequencyData

Convert phase residuals to fractional frequency via the canonical first
difference `y[k] = (x[k+1] − x[k]) / τ₀`. Length drops from `N` to
`N − 1`. The inverse direction (`_freq_to_phase`) recovers an
`N`-length record but loses the absolute phase offset, so the
phase ↔ frequency map is bijective only up to a constant: round-trip
`_freq_to_phase ∘ _phase_to_freq` returns `x[2:end] .- x[1]`.

Deviation estimators are shift-invariant (their second/third differences
annihilate the offset), so e.g. `adev(pd, m_values) ≈
adev(_phase_to_freq(pd), m_values)` up to the `N → N−1` length change.
"""
_phase_to_freq(data::PhaseData) = FrequencyData(diff(data.x) ./ data.tau0, data.tau0)

"""
    _default_m_values(N::Int, kernel::Symbol) → Vector{Int}

Octave-spaced default averaging-factor grid `[1, 2, 4, …, 2^k]` bounded
above by the kernel's algorithmic m-max — the largest `m` for which the
core L-check still yields at least one window. Used by the zero-arg
convenience methods on every public deviation API
(`adev(pd)`, `mdev(pd)`, …) so callers can skip the `m_values` argument.

Per-kernel m-max (derived from the `L`/`Ne` guard in each `_*_core`):

| `kernel`                  | m_max          |
|---------------------------|----------------|
| `:adev`, `:pdev`          | `(N − 2) ÷ 2`  |
| `:totdev`                 | `(N − 1) ÷ 2`  |
| `:mdev`, `:tdev`          | `(N − 1) ÷ 3`  |
| `:mtotdev`, `:ttotdev`    | `N ÷ 3`        |
| `:hdev`                   | `(N − 2) ÷ 3`  |
| `:htotdev`                | `(N − 1) ÷ 3`  |
| `:mhdev`, `:htdev`        | `(N − 1) ÷ 4`  |
| `:mhtotdev`               | `N ÷ 4`        |
| `:mtie`                   | `N − 1`        |

The ordinary and modified estimators require ≥2 analysis windows at the
largest τ (a single window is one raw difference, EDF ≈ 1), so their caps
sit one step below the total family, whose subsequence extension keeps a
single subsequence a valid long-τ estimate. (HTOTDEV's general branch runs
on `y = diff(x)` of length `N−1`, giving `(N−1) ÷ 3`; MTOTDEV runs on phase
directly, giving `N ÷ 3`.)

Throws `ArgumentError` for unknown kernel symbols or `N` too short to
admit any `m ≥ 1`.
"""
_default_m_values(N::Int, kernel::Symbol) = _grid(Octave, _kernel_m_max(N, kernel))
