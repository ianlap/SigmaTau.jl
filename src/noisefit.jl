# Interactive power-law fit: name a noise type and a τ range, fit its level
# with the slope fixed against what the other components leave, repeat, undo.

# The total deviations estimate the same variance as their parents.
const _NOISEFIT_RESPONSE = Dict(:adev => :adev, :totdev => :adev, :mdev => :mdev,
                                :mtotdev => :mdev, :hdev => :hdev, :htotdev => :hdev,
                                :mhdev => :mhdev, :mhtotdev => :mhdev)

"""
    NoiseFit

A power-law model built one component at a time on a measured σ(τ), by
[`noise_fit`](@ref), [`add_component!`](@ref) and [`undo!`](@ref).
[`powerlaw_sigma`](@ref), [`powerlaw_psd`](@ref) and [`clock_model`](@ref)
accept it as a model.

# Fields

$(TYPEDFIELDS)
"""
struct NoiseFit
    "The measured σ(τ)."
    result::StabilityResult
    "Estimator whose response the components are fitted with (`:adev`, `:mdev`, `:hdev` or `:mhdev`)."
    response::Symbol
    "Components in the order added: exponent α, level h_α and the τ indices it was fitted on."
    components::Vector{@NamedTuple{alpha::Int, h::Float64, range::UnitRange{Int}}}
end

Base.show(io::IO, ::MIME"text/plain", f::NoiseFit) = begin
    print(io, "NoiseFit to ", f.result.deviation_type, " (", length(f.components), " components)")
    for c in f.components
        @printf(io, "\n  h%-3s %-5s %.4g  τ %.3g–%.3g s", string(c.alpha), _nnls_symbol(c.alpha),
                c.h, f.result.tau[first(c.range)], f.result.tau[last(c.range)])
    end
end

"""
    noise_fit(r::StabilityResult) → NoiseFit

Start an interactive power-law fit on `r`, an ADEV, MDEV, HDEV or MHDEV
result or one of their total forms (which are fitted with the parent's
response). Add components with [`add_component!`](@ref) and take the last
one back with [`undo!`](@ref).

Each component's σ²(τ) is the exact response of a sampled record, as in
[`powerlaw_fit`](@ref). Its level is the weighted mean of
log(residual / response) over its own τ range, the residual being the
measured σ² less every other component, with the slope fixed and each point
weighted by its confidence interval, 1/log(hi/lo)². Without intervals the
points weigh the same. After each change every level is refitted against
the others until they settle, so the model does not depend on the order the
components were added in, and a component is not inflated by a neighbour
that was added after it.

# Examples

```julia
using SigmaTau
p = noise_gen(PhaseData, 65536, 1.0; h = Dict(2 => 1e-22, 0 => 1e-24, -2 => 1e-30))
f = noise_fit(mhdev(p))
add_component!(f, :WHPM, 1:3)      # white PM on the first three τ
add_component!(f, :WHFM, 4:9)
add_component!(f, :RWFM, 11:14)
undo!(f)                           # take random-walk FM back out
clock_model(f)
```
"""
function noise_fit(r::StabilityResult)
    resp = get(_NOISEFIT_RESPONSE, r.deviation_type, :none)
    resp === :none && throw(ArgumentError(
        "noise_fit: needs an ADEV, MDEV, HDEV or MHDEV result or a total form, got :$(r.deviation_type)"))
    isfinite(r.tau0) && r.tau0 > 0 || throw(ArgumentError("noise_fit: the result has no sample interval"))
    return NoiseFit(r, resp, @NamedTuple{alpha::Int, h::Float64, range::UnitRange{Int}}[])
end

# σ²(τₖ) of the unit component α in the fit's response.
_noisefit_unit(f::NoiseFit, α::Int, k::Int) =
    _powerlaw_var(f.response, α, round(Int, f.result.tau[k] / f.result.tau0), f.result.tau0)

"""
    noise_residual(f::NoiseFit) → Vector{Float64}

The measured σ²(τ) less the variance of every component so far. It can go
negative where the model overshoots.
"""
noise_residual(f::NoiseFit) = [f.result.dev[k]^2 - sum((c.h * _noisefit_unit(f, c.alpha, k) for c in f.components); init=0.0)
                         for k in eachindex(f.result.tau)]

# Level of power law α on `range` against the residual of every component
# but number `skip`; NaN when no τ there has a positive residual.
function _fit_level(f::NoiseFit, α::Int, range::UnitRange{Int}, skip::Int)
    r = f.result
    num = den = 0.0
    for k in range
        res = r.dev[k]^2
        for (j, c) in enumerate(f.components)
            j == skip || (res -= c.h * _noisefit_unit(f, c.alpha, k))
        end
        res > 0 && isfinite(res) || continue
        w = isempty(r.ci) ? 1.0 : 1 / log(r.ci[k].hi / r.ci[k].lo)^2
        isfinite(w) || continue
        num += w * log(res / _noisefit_unit(f, α, k))
        den += w
    end
    return den > 0 ? exp(num / den) : NaN
end

# Refit each component against the others until the levels settle. One that
# is left no positive residual keeps its level.
function _refit!(f::NoiseFit)
    for _ in 1:200
        change = 0.0
        for (j, c) in enumerate(f.components)
            h = _fit_level(f, c.alpha, c.range, j)
            isnan(h) && continue
            change = max(change, abs(log(h / c.h)))
            f.components[j] = (alpha=c.alpha, h=h, range=c.range)
        end
        change < 1e-10 && break
    end
    return f
end

"""
    add_component!(f::NoiseFit, noise, range) → NoiseFit

Add power law `noise` (an α, or a name such as `:WHFM`) to the model,
fitted at τ indices `range`, and refit the levels already there. Points
where the residual is not positive are skipped.
"""
function add_component!(f::NoiseFit, noise, range::AbstractUnitRange{<:Integer})
    α = _noise_alpha(noise)
    αmin = f.response === :adev || f.response === :mdev ? -2 : -4
    αmin <= α <= 2 || throw(ArgumentError("add_component!: α = $α does not converge in $(f.response)"))
    any(c -> c.alpha == α, f.components) && throw(ArgumentError("add_component!: α = $α is already in the model"))
    checkbounds(f.result.tau, range)
    rng = UnitRange{Int}(range)
    h = _fit_level(f, α, rng, 0)
    isnan(h) && throw(ArgumentError("add_component!: no τ in the chosen range has a positive residual"))
    push!(f.components, (alpha=α, h=h, range=rng))
    return _refit!(f)
end

"""
    undo!(f::NoiseFit) → NoiseFit

Remove the component added last and refit the rest.
"""
function undo!(f::NoiseFit)
    isempty(f.components) && throw(ArgumentError("undo!: the fit has no components"))
    pop!(f.components)
    return _refit!(f)
end

_powerlaw_tau0(f::NoiseFit) = f.result.tau0
_powerlaw_h(f::NoiseFit) = Dict(c.alpha => c.h for c in f.components)
clock_model(f::NoiseFit; kw...) = clock_model(_powerlaw_h(f); tau0=f.result.tau0, kw...)
