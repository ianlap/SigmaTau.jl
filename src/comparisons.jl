# Several clocks read against one reference: the input of every time scale
# algorithm, and a simulator that makes one with a known truth.

"""
    ClockComparisons(x::AbstractMatrix, tau0; names, reference="")
    ClockComparisons(records::AbstractVector{<:PhaseData}; names, reference="")

Phase readings of several clocks against one reference, taken at the same
instants: the input of [`timescale`](@ref).

Column ``i`` of `x` is clock ``i`` minus the reference, in seconds, with
one row per sample and `NaN` where a reading is missing. The reference may
be one of the clocks, as with a switch or a counter that measures each
against clock 1, whose column is then zero, or an outside one. `cc[i]` is
the record of clock ``i`` as [`PhaseData`](@ref), and `cc[i, j]` the
difference clock ``i`` − clock ``j``, in which the reference cancels.

# Arguments

- `x`: the readings, samples by clocks, for two or more clocks.
- `tau0`: sample interval in seconds.
- `records`: one phase record per clock, all of one length and one `tau0`.

# Keywords

- `names`: one name per clock, all different; `"clock 1"`, `"clock 2"`, …
  by default.
- `reference`: what the clocks were read against, for display.

# Fields

$(TYPEDFIELDS)

# Examples

```jldoctest
julia> using SigmaTau

julia> x = [0.0 1e-9; 0.0 3e-9; 0.0 4e-9];    # clock 2 read against clock 1

julia> cc = ClockComparisons(x, 1.0; reference = "clock 1")
ClockComparisons: 2 clocks × 3 samples at τ₀ = 1.0 s, against clock 1
  clock 1, clock 2

julia> cc[2, 1].x == x[:, 2]
true
```
"""
struct ClockComparisons
    "Phase of each clock against the reference (s), one column per clock."
    x::Matrix{Float64}
    "Sample interval τ₀ in seconds."
    tau0::Float64
    "Clock names, one per column."
    names::Vector{String}
    "What the clocks were read against."
    reference::String
    function ClockComparisons(x::AbstractMatrix{<:Real}, tau0::Real;
                              names::AbstractVector{<:AbstractString}=["clock $i" for i in 1:size(x, 2)],
                              reference::AbstractString="")
        tau0 > 0 || throw(ArgumentError("ClockComparisons: tau0 must be positive, got $tau0"))
        size(x, 2) >= 2 || throw(ArgumentError("ClockComparisons: needs two or more clocks"))
        length(names) == size(x, 2) || throw(ArgumentError("ClockComparisons: one name per clock"))
        allunique(names) || throw(ArgumentError("ClockComparisons: clock names must differ"))
        return new(Float64.(x), Float64(tau0), String.(names), String(reference))
    end
end

function ClockComparisons(records::AbstractVector{<:PhaseData};
                          names::AbstractVector{<:AbstractString}=["clock $i" for i in eachindex(records)],
                          reference::AbstractString="")
    isempty(records) && throw(ArgumentError("ClockComparisons: no records"))
    τ₀ = records[1].tau0
    all(r -> isapprox(r.tau0, τ₀; rtol=1e-9), records) ||
        throw(ArgumentError("ClockComparisons: the records must share one τ₀"))
    allequal(length(r.x) for r in records) ||
        throw(ArgumentError("ClockComparisons: the records must have the same length"))
    return ClockComparisons(reduce(hcat, (r.x for r in records)), τ₀; names, reference)
end

Base.size(c::ClockComparisons) = size(c.x)
Base.getindex(c::ClockComparisons, i::Integer) = PhaseData(c.x[:, i], c.tau0; source=c.names[i])
Base.getindex(c::ClockComparisons, i::Integer, j::Integer) =
    PhaseData(c.x[:, i] .- c.x[:, j], c.tau0; source=c.names[i] * " − " * c.names[j])

function Base.show(io::IO, ::MIME"text/plain", c::ClockComparisons)
    n, k = size(c.x)
    print(io, "ClockComparisons: ", k, " clocks × ", n, " samples at τ₀ = ", c.tau0, " s")
    isempty(c.reference) || print(io, ", against ", c.reference)
    print(io, "\n  ", join(c.names, ", "))
end

"""
    simulate_ensemble(models::AbstractVector{ClockModel}, N; reference=1, rng=nothing)
        → (ClockComparisons, Array{Float64,3})

`N` samples of an ensemble in which clock ``i`` follows `models[i]`, with
the true states behind the readings.

Each clock is drawn as by [`simulate_clock`](@ref) and read against clock
`reference`, whose own column is zero, or against a perfect outside
reference with `reference = 0`. A reading carries the white PM and the
flicker PM of its own clock's model.

# Arguments

- `models`: one [`ClockModel`](@ref) per clock, two or more, with one
  `tau0`.
- `N`: number of samples.

# Keywords

- `reference`: index of the clock the others are read against, or `0`.
- `rng`: a random number generator, such as `Xoshiro(1)`, for a
  reproducible ensemble.

# Returns

The comparisons, and the true states `truth[k, s, i]` of clock ``i`` at
sample ``k`` against true time: phase (``s = 1``), frequency (2) and, for a
model with drift, drift (3). The frequency includes the Markov states of a
model with flicker FM.

# Examples

```julia
using SigmaTau, Random
ms = [ClockModel(; tau0 = 1.0, q1 = q, q2 = 1e-6q) for q in (1e-24, 4e-24, 1.6e-23)]
cc, truth = simulate_ensemble(ms, 20_000; reference = 0, rng = Xoshiro(1))
adev(cc[2, 1])       # clock 2 − clock 1, whatever the reference
```
"""
function simulate_ensemble(models::AbstractVector{ClockModel}, N::Int; reference::Int=1, rng=nothing)
    k = length(models)
    k >= 2 || throw(ArgumentError("simulate_ensemble: needs two or more clocks"))
    0 <= reference <= k || throw(ArgumentError("simulate_ensemble: reference must be 0 or a clock index"))
    τ₀ = models[1].tau0
    all(m -> m.tau0 == τ₀, models) || throw(ArgumentError("simulate_ensemble: the models must share one τ₀"))
    truth, pm = zeros(N, 3, k), zeros(N, k)
    noise(m) = (draw = rng === nothing ? randn(N) : randn(rng, N); sqrt(_reading_var(m)) .* draw)
    for (i, m) in enumerate(models)
        _, X = simulate_clock(m, N; rng)
        np = _npoly(m)
        truth[:, 1:np, i] = X[:, 1:np]
        nf = np + length(_fm_taus(m))
        for j in np+1:size(X, 2)
            j <= nf ? (truth[:, 2, i] .+= X[:, j]) : (pm[:, i] .+= X[:, j])
        end
    end
    ref = reference == 0 ? zeros(N) : truth[:, 1, reference]
    x = reduce(hcat, (i == reference ? zeros(N) : truth[:, 1, i] .- ref .+ pm[:, i] .+ noise(models[i]) for i in 1:k))
    names = ["clock $i" for i in 1:k]
    return ClockComparisons(x, τ₀; names, reference=reference == 0 ? "ideal reference" : names[reference]), truth
end

"""
    clock_models(cc::ClockComparisons; alphas=[2, 1, 0, -2], flicker=false)
        → Vector{ClockModel}

One [`ClockModel`](@ref) per clock, fitted from the comparisons alone, for
the Kalman time scales of [`timescale`](@ref).

Each clock's Allan deviation is separated from the others' by the
Groslambert covariance of the clock differences ([`groslambert`](@ref)),
in which the reference cancels. [`powerlaw_fit`](@ref) fits the components
`alphas` to that curve, and [`clock_model`](@ref) turns the fit into a
model. This takes three or more clocks.

The default components are white PM, flicker PM, white FM and random-walk
FM, and the flicker PM is folded into `R`. Flicker PM needs its own term
in the fit: without one it passes for white FM and inflates `q1`. With
`flicker = true`, flicker FM joins the fit and each model keeps its
flicker FM and flicker PM as Markov states. At the longest τ of a short
record the levels are poorly separated, so compare each model with its
curve. A clock whose curve leaves the fit too few positive points throws an
error.

# Arguments

- `cc`: the comparisons, of three or more clocks.

# Keywords

- `alphas`: exponents ``α`` of the components fitted.
- `flicker`: fit flicker FM as well and keep both flickers as Markov
  states.

# Examples

```julia
using SigmaTau, Random
ms = [ClockModel(; tau0 = 1.0, q1 = q, q2 = 1e-6q) for q in (1e-24, 4e-24, 1.6e-23)]
cc, _ = simulate_ensemble(ms, 20_000; reference = 0, rng = Xoshiro(1))
fitted = clock_models(cc)
[m.q1 for m in fitted]        # near 1e-24, 4e-24, 1.6e-23
```
"""
function clock_models(cc::ClockComparisons; alphas=[2, 1, 0, -2], flicker::Bool=false)
    k = size(cc.x, 2)
    k >= 3 || throw(ArgumentError("clock_models: needs three or more clocks, got $k"))
    pair = Matrix{PhaseData{Float64}}(undef, k, k)
    for i in 1:k, j in (i + 1):k
        pair[i, j] = cc[i, j]
    end
    flicker && (alphas = union(alphas, [1, -1]))
    return ClockModel[clock_model(powerlaw_fit(r; alphas); flicker) for r in groslambert(pair)]
end
