# Averaging-factor grids (TauMode, tau_values) and phase/frequency conversion.
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
and `_default_m_values`.

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
| `:mtie`, `:tierms`        | `N − 1`        |
| `:theo1`                  | `2·((N − 1) ÷ 2)` (largest even ≤ N − 1) |
| `:theoh`                  | `⌊0.75 · 2·((N − 1) ÷ 2)⌋` (τ-grid units) |

The ordinary and modified estimators require ≥2 analysis windows at the largest
τ (a single window is one raw difference, EDF ≈ 1), so their caps sit one step
below the total family, whose subsequence extension keeps a single subsequence a
valid long-τ estimate. (HTOTDEV's general branch runs on `y = diff(x)` of length
`N−1`, giving `(N−1) ÷ 3`; MTOTDEV runs on phase directly, giving `N ÷ 3`.)
Thêo1 is defined at even `m ≤ N − 1` (SP1065 eq. 30); ThêoH's grid is in τ
units (`τ = m·τ0`), reaching Thêo1's effective-τ limit `0.75·(N−1)·τ0`.

Throws `ArgumentError` for an unknown kernel symbol or an `N` too short to admit
any `m ≥ 1`.
"""
function _kernel_m_max(N::Int, kernel::Symbol)
    # Ordinary/modified estimators use raw differences, so the largest τ must
    # retain ≥2 analysis windows; a single window is one difference (EDF ≈ 1)
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
    elseif kernel === :mtie || kernel === :tierms
        N - 1
    elseif kernel === :theo1
        2 * ((N - 1) ÷ 2)
    elseif kernel === :theoh
        floor(Int, 0.75 * (2 * ((N - 1) ÷ 2)))
    else
        throw(ArgumentError("tau_values: unknown kernel symbol :$kernel"))
    end
    m_max < 1 && throw(ArgumentError(
        "tau_values: N=$N is too short to support any m for :$kernel"))
    return m_max
end

"""
    _grid(mode::TauMode, m_max::Int) → Vector{Int}

Averaging-factor grid for `mode`, clamped to `[1, m_max]`. `Octave` uses exact
powers of two; the other geometric modes round their factor powers to `Int`,
deduplicate, and clamp.
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
(`:adev`, `:mdev`, …) that determines the m-max clamp; see `_kernel_m_max`.

For `:theo1` (defined at even averaging factors only, SP1065 eq. 30) odd grid
entries are rounded up to the next even value, deduplicated, and clamped to the
even m-max. `:theoh` grids are plain integers: its m-values are τ-grid factors
(target `τ = m·τ0`) and the even-factor conversion for the Thêo1 segment
happens inside [`theoh`](@ref).

`tau_values(Octave, N, kernel)` equals `_default_m_values``(N, kernel)`.
"""
function tau_values(mode::TauMode, N::Int, kernel::Symbol)
    m_max = _kernel_m_max(N, kernel)
    ms = _grid(mode, m_max)
    kernel === :theo1 && return unique!([min(2 * ((m + 1) ÷ 2), m_max) for m in ms])
    return ms
end


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
integral `x[1] = 0`, `x[k+1] = τ₀ · Σⱼ₌₁ᵏ y[j]` (NIST SP1065 §1.6; the
conversion Stable32 and allantools perform). `N` frequency samples give
`N + 1` phase samples, so every `y[k]` survives as the phase difference
`x[k+1] − x[k]` and `adev(FrequencyData(y, τ₀), …)` equals
`adev(PhaseData([0; cumsum(y)]·τ₀, τ₀), …)`. This is the exact inverse of
the first difference up to the phase offset `x[1]`.
"""
function _freq_to_phase(data::FrequencyData)
    y = data.y
    x = Vector{Float64}(undef, length(y) + 1)
    x[1] = 0.0
    acc = 0.0
    @inbounds for k in eachindex(y)
        acc += y[k]
        x[k+1] = acc * data.tau0
    end
    return PhaseData(x, data.tau0)
end

"""
    FrequencyData(p::PhaseData) → FrequencyData
    PhaseData(f::FrequencyData) → PhaseData

The record as the other kind, as Stable32's Convert function and allantools'
`phase2frequency` / `frequency2phase` do [riley-2008-sp1065](@cite). Phase
becomes frequency by the first difference `y[k] = (x[k+1] − x[k]) / τ₀`, so
`N` phase samples give `N − 1` frequency samples and a missing phase sample
makes the two beside it missing. Frequency becomes phase by the running sum
`x = τ₀ · [0; cumsum(y)]`, `N + 1` samples from 0; a missing frequency sample
leaves every later phase unknown, so fill it first with [`fillgaps`](@ref).
`τ₀`, `source` and `start` carry over.

The estimators take either kind and convert internally, so this is for
plotting, saving or handing the record to another tool. The round trip
`PhaseData(FrequencyData(p))` returns `p.x .- p.x[1]`.

```jldoctest
julia> using SigmaTau

julia> FrequencyData(PhaseData([0.0, 1.0, 3.0], 0.5)).y
2-element Vector{Float64}:
 2.0
 4.0
```
"""
FrequencyData(p::PhaseData) = FrequencyData(diff(p.x) ./ p.tau0, p.tau0; source=p.source, start=p.start)

function PhaseData(f::FrequencyData)
    any(isnan, f.y) && throw(ArgumentError(
        "PhaseData: a missing frequency sample leaves the phase after it unknown; fill the gaps first (fillgaps)"))
    return PhaseData(_freq_to_phase(f).x, f.tau0; source=f.source, start=f.start)
end

# Either record as the other kind, and its raw samples.
_as_phase(d::PhaseData) = d
_as_phase(d::FrequencyData) = _freq_to_phase(d)
_as_freq(d::FrequencyData) = d
_as_freq(d::PhaseData) = FrequencyData(d)
_samples(d::PhaseData) = d.x
_samples(d::FrequencyData) = d.y

"""
    _default_m_values(N::Int, kernel::Symbol) → Vector{Int}

Octave-spaced averaging-factor grid `[1, 2, 4, …]` capped at the kernel's
algorithmic m-max (see `_kernel_m_max`); what every zero-argument
deviation call (`adev(pd)`, `mdev(pd)`, …) uses. Equivalent to
`tau_values(Octave, N, kernel)`.
"""
_default_m_values(N::Int, kernel::Symbol) = tau_values(Octave, N, kernel)
