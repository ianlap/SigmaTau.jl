# Several clocks read against one reference: the input of every timescale
# algorithm, and a simulator that makes one with a known truth.

"""
    ClockComparisons(x::Matrix, tau0; names, reference="")
    ClockComparisons(records::Vector{<:PhaseData}; names, reference="")

Phase readings of several clocks against one reference, sampled together:
column i of `x` is clock i minus the reference, in seconds, one row per
sample, `NaN` where a reading is missing. The reference may be one of the
clocks (a switch or counter measuring each against clock 1, whose column is
then zero) or an outside one, named by `reference`.

`comparisons[i]` is clock i's record as `PhaseData`, and
`comparisons[i, j]` the difference clock i − clock j, which no longer
depends on the reference.

$(TYPEDFIELDS)
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
    simulate_ensemble(models::Vector{ClockModel}, N; reference=1, rng=nothing)
        → (ClockComparisons, Array{Float64,3})

`N` samples of an ensemble whose clock i follows `models[i]`, read against
clock `reference` (or, with `reference = 0`, against a perfect outside
reference). Each reading carries its clock's white PM `R`. Returns the
comparisons and the true states, `truth[k, s, i]` for sample k, state s
and clock i (phase, frequency, and drift for a model with drift; zero
beyond a clock's own states). The models must share one τ₀.
"""
function simulate_ensemble(models::AbstractVector{ClockModel}, N::Int; reference::Int=1, rng=nothing)
    k = length(models)
    k >= 2 || throw(ArgumentError("simulate_ensemble: needs two or more clocks"))
    0 <= reference <= k || throw(ArgumentError("simulate_ensemble: reference must be 0 or a clock index"))
    τ₀ = models[1].tau0
    all(m -> m.tau0 == τ₀, models) || throw(ArgumentError("simulate_ensemble: the models must share one τ₀"))
    truth = zeros(N, 3, k)
    noise(m) = (draw = rng === nothing ? randn(N) : randn(rng, N); sqrt(m.R) .* draw)
    for (i, m) in enumerate(models)
        _, X = simulate_clock(ClockModel(m.tau0, 0.0, m.q1, m.q2, m.q3, m.drift), N; rng)
        truth[:, 1:size(X, 2), i] = X
    end
    ref = reference == 0 ? zeros(N) : truth[:, 1, reference]
    x = reduce(hcat, (i == reference ? zeros(N) : truth[:, 1, i] .- ref .+ noise(models[i]) for i in 1:k))
    names = ["clock $i" for i in 1:k]
    return ClockComparisons(x, τ₀; names, reference=reference == 0 ? "ideal reference" : names[reference]), truth
end

"""
    clock_models(cc::ClockComparisons; alphas=[2, 0, -2]) → Vector{ClockModel}

One [`ClockModel`](@ref) per clock, fitted from the comparisons alone, for the
Kalman timescales (`:greenhall`, `:jones_tryon`), which need one. Needs three
or more clocks: each clock's Allan deviation comes from the Groslambert
covariance of the clock differences ([`groslambert`](@ref)), so the reference
and the other clocks' noise drop out, and the power-law fit of
[`powerlaw_fit`](@ref) on that curve, with the components `alphas` (white PM,
white FM and random-walk FM by default; the Allan deviation cannot tell flicker
PM from white PM), gives the model. A clock whose curve the fit cannot use
(too few positive points) throws.
"""
function clock_models(cc::ClockComparisons; alphas=[2, 0, -2])
    k = size(cc.x, 2)
    k >= 3 || throw(ArgumentError("clock_models: needs three or more clocks, got $k"))
    pair = Matrix{PhaseData{Float64}}(undef, k, k)
    for i in 1:k, j in (i + 1):k
        pair[i, j] = cc[i, j]
    end
    return ClockModel[clock_model(powerlaw_fit(r; alphas)) for r in groslambert(pair)]
end
