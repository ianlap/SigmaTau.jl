# Averaging-factor grids (TauMode, tau_values) and phase/frequency conversion.
"""
    TauMode

Spacing of the averaging-factor grid. Pass an instance in place of
`m_values`, as in `adev(data, Decade)`, or to [`stability`](@ref) and
[`tau_values`](@ref).

| Instance        | Ratio of successive factors | Grid                            |
|:----------------|:----------------------------|:--------------------------------|
| `AllTaus`       | every integer               | `1, 2, 3, 4, 5, …`              |
| `Octave`        | 2                           | `1, 2, 4, 8, 16, …`             |
| `HalfOctave`    | √2                          | `1, 2, 3, 4, 6, 8, 11, 16, …`   |
| `QuarterOctave` | 2^(1/4)                     | `1, 2, 3, 4, 5, 6, 7, 8, 10, …` |
| `Decade`        | 10                          | `1, 10, 100, 1000, …`           |
| `HalfDecade`    | √10                         | `1, 3, 10, 32, 100, …`          |

`Octave` is the default of every estimator. The geometric grids are rounded
to integers, duplicates are dropped, and each grid stops at the largest
factor the estimator supports for the record length.
"""
@enum TauMode AllTaus Octave HalfOctave QuarterOctave Decade HalfDecade

"""
    _kernel_m_max(N::Int, kernel::Symbol) → Int

Largest averaging factor `m` that `kernel` supports on `N` phase samples.

| `kernel`               | Largest `m`                                     |
|:-----------------------|:------------------------------------------------|
| `:adev`, `:pdev`       | `(N − 2) ÷ 2`                                   |
| `:totdev`              | `(N − 1) ÷ 2`                                   |
| `:mdev`, `:tdev`       | `(N − 1) ÷ 3`                                   |
| `:mtotdev`, `:ttotdev` | `N ÷ 3`                                         |
| `:hdev`                | `(N − 2) ÷ 3`                                   |
| `:htotdev`             | `(N − 1) ÷ 3`                                   |
| `:mhdev`, `:htdev`     | `(N − 1) ÷ 4`                                   |
| `:mhtotdev`            | `N ÷ 4`                                         |
| `:mtie`, `:tierms`     | `N − 1`                                         |
| `:theo1`               | `2·((N − 1) ÷ 2)`, the largest even `m ≤ N − 1` |
| `:theoh`               | `⌊0.75 · 2·((N − 1) ÷ 2)⌋`, in τ-grid units      |

The Allan, Hadamard and parabolic estimators need two terms at the largest
τ, since one term is a single difference with an EDF near 1. The total
estimators keep a single subsequence, which its extension makes a usable
estimate. HTOTDEV works on the `N − 1` frequency values, hence
`(N − 1) ÷ 3` against MTOTDEV's `N ÷ 3`. ThêoH's grid is in units of
`τ = m·τ0` and reaches Thêo1's largest averaging time, `0.75·(N − 1)·τ0`.

Throws an `ArgumentError` for an unknown kernel or an `N` too short for
`m = 1`.
"""
function _kernel_m_max(N::Int, kernel::Symbol)
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

Averaging factors of `mode` up to `m_max`. `Octave` uses exact powers of
two; the other geometric modes round the powers of their ratio and drop
duplicates.
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

Averaging factors ``m`` of the grid `mode` for a record of `N` phase
samples, up to the largest that the estimator `kernel` (`:adev`, `:mdev`,
…) supports. Despite the name, the values are averaging factors, not times:
``τ = mτ_0``. An estimator called with a [`TauMode`](@ref) uses this grid,
so `adev(data, Decade)` is `adev(data, tau_values(Decade, N, :adev))`.

For `:theo1`, which takes even factors only, odd entries are rounded up to
the next even value. A `:theoh` grid is in units of ``τ = mτ_0``, and
[`theoh`](@ref) chooses the Thêo1 factor for each point.

# Examples

```jldoctest
julia> using SigmaTau

julia> print(tau_values(Decade, 1000, :adev))
[1, 10, 100]

julia> print(tau_values(Octave, 20, :theo1))
[2, 4, 8, 16]
```
"""
function tau_values(mode::TauMode, N::Int, kernel::Symbol)
    m_max = _kernel_m_max(N, kernel)
    ms = _grid(mode, m_max)
    kernel === :theo1 && return unique!([min(2 * ((m + 1) ÷ 2), m_max) for m in ms])
    return ms
end


"""
    _f64(v::AbstractVector{<:Real}) → Vector{Float64}

The samples as the `Vector{Float64}` that the kernels and `identify_noise`
take: `v` itself when it already is one, otherwise a converted copy. This
lets a record of any real element type, such as `PhaseData{Float32}`, pass
through the estimators.
"""
_f64(v::Vector{Float64}) = v
_f64(v::AbstractVector{<:Real}) = Vector{Float64}(v)

"""
    _freq_to_phase(data::FrequencyData) → PhaseData

Phase from fractional frequency by the running sum `x[1] = 0`,
`x[k+1] = x[k] + τ₀·y[k]` (SP1065 §3.4; the Stable32 Convert function).
`N` frequency samples give `N + 1` phase samples, so each `y[k]` survives
as the difference `(x[k+1] − x[k])/τ₀`, and the conversion inverts the
first difference up to the phase offset `x[1]`.
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

The record as the other kind
([SP1065 §3.4](@cite riley-2008-sp1065); the Convert function of the
[Stable32 User Manual](@cite hamilton-2008-stable32-manual)).

Phase becomes frequency by the first difference
``y_k = (x_{k+1} - x_k)/τ_0``: ``N`` phase samples give ``N - 1`` frequency
samples, and a missing phase sample makes the two frequency samples beside
it missing. Frequency becomes phase by the running sum
``x_{k+1} = x_k + τ_0 y_k`` from ``x_1 = 0``, which gives ``N + 1`` samples.
A missing frequency sample leaves every later phase unknown, so fill it
first with [`fillgaps`](@ref). `tau0`, `source` and `start` carry over.

The estimators take either kind and convert internally, so the conversion
is for plotting, saving or handing a record to another tool. The round trip
`PhaseData(FrequencyData(p))` returns `p.x .- p.x[1]`. The mean frequency
is kept: allantools' `frequency2phase` removes it before integrating, which
changes MTIE and TIE rms and no other estimator.

# Examples

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

The grid of an estimator called without `m_values`:
`tau_values(Octave, N, kernel)`.
"""
_default_m_values(N::Int, kernel::Symbol) = tau_values(Octave, N, kernel)
