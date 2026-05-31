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

| `kernel`                                     | m_max          |
|----------------------------------------------|----------------|
| `:adev`, `:totdev`, `:pdev`                  | `(N − 1) ÷ 2`  |
| `:mdev`, `:tdev`, `:mtotdev`, `:ttotdev`     | `N ÷ 3`        |
| `:hdev`, `:htotdev`                          | `(N − 1) ÷ 3`  |
| `:mhdev`, `:htdev`, `:mhtotdev`              | `N ÷ 4`        |
| `:mtie`                                      | `N − 1`        |

(HTOTDEV's general branch operates on `y = diff(x)` of length `N−1`, so its
`n_iter = (N−1) − 3m + 1 ≥ 1` constraint matches HDEV's even though MTOTDEV —
which runs on phase directly — uses `N ÷ 3`.)

Throws `ArgumentError` for an unknown kernel symbol or an `N` too short to admit
any `m ≥ 1`.
"""
function _kernel_m_max(N::Int, kernel::Symbol)
    m_max = if kernel === :adev || kernel === :totdev || kernel === :pdev
        (N - 1) ÷ 2
    elseif kernel === :mdev || kernel === :tdev ||
           kernel === :mtotdev || kernel === :ttotdev
        N ÷ 3
    elseif kernel === :hdev || kernel === :htotdev
        (N - 1) ÷ 3
    elseif kernel === :mhdev || kernel === :htdev || kernel === :mhtotdev
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
