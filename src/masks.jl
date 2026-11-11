# ITU-T wander masks for MTIE and TDEV: the wander generation limits of
# G.811 (primary reference clock), G.813 (SDH equipment clock) and G.8262
# (synchronous equipment clock), and the G.823 network limits at
# synchronization interfaces of the 2048 kbit/s hierarchy.

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
    TELECOM_MASKS

Names of the ITU-T wander masks that [`telecom_mask`](@ref) and
[`mask_check`](@ref) know. Each limits MTIE and TDEV.

- `:G811`: wander generation of a primary reference clock
  ([G.811, clause 6.1](@cite ITU-G811)).
- `:G813_opt1`, `:G813_opt2`: wander generation of an SDH equipment clock
  in locked mode at constant temperature
  ([G.813, clause 7.1](@cite ITU-G813)). Option 1 is for networks of the
  2048 kbit/s hierarchy, option 2 for the 1544 kbit/s hierarchy.
- `:G8262_opt1`, `:G8262_opt2`: the same for a synchronous equipment
  clock, with the limits of G.813
  ([G.8262, clause 8.1](@cite ITU-G8262)).
- `:G823_PRC`, `:G823_SSU`, `:G823_SEC`, `:G823_PDH`: network limits for
  wander at PRC, SSU, SEC and PDH synchronization interfaces
  ([G.823, clause 6.2](@cite ITU-G823)).

The clock recommendations measure MTIE and TDEV through a 10 Hz
first-order low-pass filter at a sample interval of 1/30 s or less, and
TDEV over a record of at least ``12τ``. This package applies no such
filter: the record is compared as it was measured.
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
    telecom_mask(mask, stat, τ) → Float64 or Vector{Float64}

Limit of `mask`, one of [`TELECOM_MASKS`](@ref), for `stat`, `:mtie` or
`:tdev`, at the observation interval `τ` in seconds, or at each of a vector
of them. The limit is in seconds, and `NaN` outside the range of ``τ`` the
recommendation covers.

# Examples

```jldoctest
julia> using SigmaTau

julia> round(telecom_mask(:G8262_opt1, :tdev, 10.0) * 1e9; digits=6)   # in ns
3.2
```
"""
telecom_mask(mask::Symbol, stat::Symbol, τ::Real) = _mask_limit(_mask_pieces(mask, stat), Float64(τ))
telecom_mask(mask::Symbol, stat::Symbol, τ::AbstractVector{<:Real}) =
    (p = _mask_pieces(mask, stat); Float64[_mask_limit(p, Float64(t)) for t in τ])

"""
    mask_check(r::StabilityResult, mask) → NamedTuple

Compare an [`mtie`](@ref) or [`tdev`](@ref) result with a telecom mask, one
of [`TELECOM_MASKS`](@ref). The result has the fields

- `mask`: the name of the mask.
- `tau`: the observation intervals of `r`.
- `limit`: the mask at each of them; `NaN` where the mask does not reach.
- `pass`: `true` when every value is at or under its limit.
- `margin`: the smallest ratio of limit to value. A margin above 1 passes,
  and the reciprocal of a smaller one is the factor by which the record
  exceeds the mask.

Points outside the range of the mask, and points without a finite value,
are not compared. With none left, `margin` is `NaN` and `pass` is `false`.

# Examples

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
