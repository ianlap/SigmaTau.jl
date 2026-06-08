"""
    StabilitySuite

Ordered collection of [`StabilityResult`](@ref)s computed from one timing record
in a single [`stability`](@ref) call, plus the session metadata needed to
reproduce and serialize the analysis.

Index by deviation symbol to recover a result (`suite[:adev]`) or by position
(`suite[1]`). The suite is iterable in request order and supports `keys`,
`haskey`, `length`, and `pairs`, so a CLI or GUI can render it as an ordered
table.

$(TYPEDFIELDS)
"""
struct StabilitySuite
    "Per-deviation results, in the order they were requested."
    results::Vector{StabilityResult}
    "Input record kind: `:phase` or `:frequency`."
    data_kind::Symbol
    "Base sample interval τ₀ in seconds."
    tau0::Float64
    "Number of input samples N."
    n::Int
    "Confidence level used for the CI bounds (`NaN` when `ci=false`)."
    confidence::Float64
    "Tau-grid selector, recorded as a Symbol (e.g. `:Octave`, or `:explicit`)."
    tau_mode::Symbol
end

Base.length(s::StabilitySuite) = length(s.results)
Base.iterate(s::StabilitySuite, state...) = iterate(s.results, state...)
Base.eltype(::Type{StabilitySuite}) = StabilityResult
Base.keys(s::StabilitySuite) = Symbol[r.deviation_type for r in s.results]
Base.haskey(s::StabilitySuite, k::Symbol) = any(r -> r.deviation_type === k, s.results)
Base.pairs(s::StabilitySuite) = (r.deviation_type => r for r in s.results)

Base.getindex(s::StabilitySuite, i::Integer) = s.results[i]

function Base.getindex(s::StabilitySuite, k::Symbol)
    idx = findfirst(r -> r.deviation_type === k, s.results)
    idx === nothing && throw(KeyError(k))
    return s.results[idx]
end

function Base.show(io::IO, s::StabilitySuite)
    devs = join(("$d" for d in keys(s)), ", ")
    ci = isnan(s.confidence) ? "no CI" : "CI@$(s.confidence)"
    print(io, "StabilitySuite(", length(s), " devs [", devs, "], ",
              s.data_kind, ", N=", s.n, ", τ₀=", s.tau0, " s, ", ci, ")")
end

function Base.show(io::IO, ::MIME"text/plain", s::StabilitySuite)
    ci = isnan(s.confidence) ? "no CI" : "CI @ $(s.confidence)"
    print(io, "StabilitySuite: ", length(s), " deviation(s), ", s.data_kind,
              " data, N=", s.n, ", τ₀=", s.tau0, " s, taus=", s.tau_mode, ", ", ci)
    for r in s.results
        n = length(r.tau)
        rng = n == 0 ? "" : "  τ∈[$(r.tau[1]), $(r.tau[end])] s"
        print(io, "\n  ", rpad(string(r.deviation_type), 9), lpad(n, 3), " pts", rng)
    end
end
