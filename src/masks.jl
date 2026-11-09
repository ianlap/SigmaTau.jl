# ITU-T wander masks for MTIE and TDEV: the clock wander-generation limits of
# G.811 (PRC), G.813 (SDH equipment clock) and G.8262 (synchronous Ethernet
# equipment clock), and the G.823 network limits at 2048 kbit/s
# synchronization interfaces.

# One piece of a mask: for lo < τ ≤ hi, limit = a + b·τ^p + c·τ in ns.
const _MaskPiece = NTuple{6,Float64}

_pieces(rows::Tuple...) = _MaskPiece[(Float64(r[1]), Float64(r[2]), Float64(r[3]), Float64(r[4]),
                                      Float64(r[5]), Float64(r[6])) for r in rows]

# Each mask: a description, then the MTIE and TDEV pieces (lo, hi, a, b, p, c).
const _G811_MTIE = _pieces((0.1, 1000, 25, 0, 0, 0.275), (1000, Inf, 290, 0, 0, 0.01))
const _G811_TDEV = _pieces((0.1, 100, 3, 0, 0, 0), (100, 1000, 0, 0, 0, 0.03), (1000, 10_000, 30, 0, 0, 0))
const _EC1_MTIE = _pieces((0.1, 1, 40, 0, 0, 0), (1, 100, 0, 40, 0.1, 0), (100, 1000, 0, 25.25, 0.2, 0))
const _EC1_TDEV = _pieces((0.1, 25, 3.2, 0, 0, 0), (25, 100, 0, 0.64, 0.5, 0), (100, 1000, 6.4, 0, 0, 0))
const _EC2_MTIE = _pieces((0.1, 1, 20, 0, 0, 0), (1, 10, 0, 20, 0.48, 0), (10, 1000, 60, 0, 0, 0))
const _EC2_TDEV = _pieces((0.1, 2.5, 0, 3.2, -0.5, 0), (2.5, 40, 2, 0, 0, 0), (40, 1000, 0, 0.32, 0.5, 0),
                          (1000, 10_000, 10, 0, 0, 0))
# G.823's long-τ terms: 433·τ^0.2 + 0.01·τ (MTIE), 58 + 1.2·τ^0.5 + 0.0003·τ (TDEV).
_g823_mtie(t1, m1, c, t2) = _pieces((0.1, t1, m1, 0, 0, 0), (t1, t2, 0, 0, 0, c),
                                 (t2, 2000, 2000, 0, 0, 0), (2000, Inf, 0, 433, 0.2, 0.01))
_g823_tdev(t1, d1) = _pieces((0.1, t1, d1, 0, 0, 0), (t1, 100, 0, 0, 0, 0.7), (100, 1e6, 58, 1.2, 0.5, 0.0003))

const _MASKS = Dict{Symbol,Tuple{String,Vector{_MaskPiece},Vector{_MaskPiece}}}(
    :G811       => ("ITU-T G.811 primary reference clock", _G811_MTIE, _G811_TDEV),
    :G813_opt1  => ("ITU-T G.813 SDH equipment clock, option 1", _EC1_MTIE, _EC1_TDEV),
    :G813_opt2  => ("ITU-T G.813 SDH equipment clock, option 2", _EC2_MTIE, _EC2_TDEV),
    :G8262_opt1 => ("ITU-T G.8262 synchronous Ethernet equipment clock, option 1", _EC1_MTIE, _EC1_TDEV),
    :G8262_opt2 => ("ITU-T G.8262 synchronous Ethernet equipment clock, option 2", _EC2_MTIE, _EC2_TDEV),
    :G823_PRC   => ("ITU-T G.823 network limit at PRC interfaces", _G811_MTIE,
                    [_G811_TDEV; _pieces((10_000, 1e6, 27, 0, 0, 0.0003))]),
    :G823_SSU   => ("ITU-T G.823 network limit at SSU interfaces", _g823_mtie(2.5, 25, 10, 200), _g823_tdev(4.3, 3)),
    :G823_SEC   => ("ITU-T G.823 network limit at SEC interfaces", _g823_mtie(2.5, 250, 100, 20), _g823_tdev(17.14, 12)),
    :G823_PDH   => ("ITU-T G.823 network limit at PDH synchronization interfaces", _g823_mtie(7.3, 732, 100, 20),
                    _g823_tdev(48, 34)),
)

"""
Names of the telecom wander masks [`telecom_mask`](@ref) and
[`mask_check`](@ref) know, each with its MTIE and TDEV limit:

- `:G811`: primary reference clock, ITU-T G.811 [ITU-G811](@cite).
- `:G813_opt1`, `:G813_opt2`: SDH equipment clock (SEC), ITU-T G.813
  [ITU-G813](@cite), options 1 (2048 kbit/s networks) and 2 (1544 kbit/s).
- `:G8262_opt1`, `:G8262_opt2`: synchronous Ethernet equipment clock
  (EEC), ITU-T G.8262 [ITU-G8262](@cite); the same limits as G.813.
- `:G823_PRC`, `:G823_SSU`, `:G823_SEC`, `:G823_PDH`: ITU-T G.823 network
  limits for wander at PRC, SSU, SEC and PDH synchronization interfaces
  [ITU-G823](@cite).

The clock masks are wander generation in the locked mode at constant
temperature; the recommendations measure through a 10 Hz first-order
low-pass filter, sampling at 1/30 s or faster.
"""
const TELECOM_MASKS = (:G811, :G813_opt1, :G813_opt2, :G8262_opt1, :G8262_opt2,
                       :G823_PRC, :G823_SSU, :G823_SEC, :G823_PDH)

function _mask_pieces(mask::Symbol, stat::Symbol)
    m = get(_MASKS, mask, nothing)
    m === nothing && throw(ArgumentError("unknown telecom mask :" * String(mask) *
                                         "; use one of " * join(map(String, TELECOM_MASKS), ", ")))
    stat === :mtie && return m[2]
    stat === :tdev && return m[3]
    throw(ArgumentError("telecom masks limit :mtie and :tdev, not :" * String(stat)))
end

function _mask_limit(pieces::Vector{_MaskPiece}, τ::Float64)
    for (lo, hi, a, b, p, c) in pieces
        lo < τ <= hi && return (a + b * τ^p + c * τ) * 1e-9
    end
    return NaN
end

"""
    telecom_mask(mask, stat, τ) → limit in seconds

The limit of `mask` (one of [`TELECOM_MASKS`](@ref)) for `stat` (`:mtie`
or `:tdev`) at observation interval `τ` in seconds, or a vector of limits for
a vector of `τ`. `NaN` outside the τ range the recommendation covers.

```jldoctest
julia> round(telecom_mask(:G8262_opt1, :tdev, 10.0) * 1e9; digits=6)   # ns
3.2
```
"""
telecom_mask(mask::Symbol, stat::Symbol, τ::Real) = _mask_limit(_mask_pieces(mask, stat), Float64(τ))
telecom_mask(mask::Symbol, stat::Symbol, τ::AbstractVector{<:Real}) =
    (p = _mask_pieces(mask, stat); Float64[_mask_limit(p, Float64(t)) for t in τ])

"""
    mask_check(r::StabilityResult, mask) → (; mask, tau, limit, pass, margin)

Compare an [`mtie`](@ref) or [`tdev`](@ref) result with a telecom mask (one
of [`TELECOM_MASKS`](@ref)). `limit[i]` is the mask at `r.tau[i]` (`NaN`
where the mask does not reach). `pass` is true when every value is at or
under its limit, and `margin` is the smallest `limit / value`: above 1
passes, and its reciprocal is how far a failing record is over. Points
outside the mask's τ range, or without a finite value, are not compared;
with none left, `margin` is `NaN` and `pass` is false.

```julia
r = mtie(read_phase("eec.txt"))
c = mask_check(r, :G8262_opt1)
c.pass, c.margin
```
"""
function mask_check(r::StabilityResult, mask::Symbol)
    r.deviation_type in (:mtie, :tdev) ||
        throw(ArgumentError("mask_check: telecom masks limit MTIE and TDEV, not :" * String(r.deviation_type)))
    limit = telecom_mask(mask, r.deviation_type, r.tau)
    margin = NaN
    for i in eachindex(limit)
        (isfinite(limit[i]) && isfinite(r.dev[i])) || continue
        q = r.dev[i] > 0 ? limit[i] / r.dev[i] : Inf
        margin = isnan(margin) ? q : min(margin, q)
    end
    return (mask=mask, tau=r.tau, limit=limit, pass=margin >= 1, margin=margin)
end
