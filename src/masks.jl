# ITU-T wander masks for MTIE and TDEV: the wander generation limits of the
# clocks of G.811, G.811.1, G.812, G.813, G.8262, G.8262.1, G.8272 and
# G.8272.1, and the G.823 network limits at synchronization interfaces of
# the 2048 kbit/s hierarchy. A further mask is one more row of _MASK_TABLE.

# One piece of a mask: for lo < τ ≤ hi, limit = a + b·τ^p + c·τ in ns.
const _MaskPiece = NTuple{6,Float64}

_pieces(rows::Tuple...) = _MaskPiece[(Float64(r[1]), Float64(r[2]), Float64(r[3]), Float64(r[4]),
                                      Float64(r[5]), Float64(r[6])) for r in rows]

# The bound just under x: as `hi` it leaves x out (τ < x), as `lo` it takes x in (x ≤ τ).
_under(x) = prevfloat(Float64(x))

# A mask: what it limits, the edition and tables its numbers were read from,
# and the MTIE and TDEV pieces (lo, hi, a, b, p, c).
struct _Mask
    title::String
    source::String
    mtie::Vector{_MaskPiece}
    tdev::Vector{_MaskPiece}
end

const _G811_MTIE = _pieces((0.1, 1000, 25, 0, 0, 0.275), (1000, Inf, 290, 0, 0, 0.01))
const _EC1_MTIE = _pieces((0.1, 1, 40, 0, 0, 0), (1, 100, 0, 40, 0.1, 0), (100, 1000, 0, 25.25, 0.2, 0))
const _EC1_TDEV = _pieces((0.1, 25, 3.2, 0, 0, 0), (25, 100, 0, 0.64, 0.5, 0), (100, 1000, 6.4, 0, 0, 0))
const _EC2_MTIE = _pieces((0.1, 1, 20, 0, 0, 0), (1, 10, 0, 20, 0.48, 0), (10, 1000, 60, 0, 0, 0))
const _EC2_TDEV = _pieces((0.1, 2.5, 0, 3.2, -0.5, 0), (2.5, 40, 2, 0, 0, 0), (40, 1000, 0, 0.32, 0.5, 0),
                          (1000, 10_000, 10, 0, 0, 0))
# G.823's long-τ terms: 433·τ^0.2 + 0.01·τ (MTIE), 58 + 1.2·τ^0.5 + 0.0003·τ (TDEV).
_g823_mtie(t1, m1, c, t2) = _pieces((0.1, t1, m1, 0, 0, 0), (t1, t2, 0, 0, 0, c),
                                 (t2, 2000, 2000, 0, 0, 0), (2000, Inf, 0, 433, 0.2, 0.01))
_g823_tdev(t1, d1) = _pieces((0.1, t1, d1, 0, 0, 0), (t1, 100, 0, 0, 0, 0.7), (100, 1e6, 58, 1.2, 0.5, 0.0003))
# PRTC-A and PRTC-B share the slope of G.811 and differ in the ceiling it meets at `t`.
_prtc_mtie(t, ceiling) = _pieces((0.1, t, 25, 0, 0, 0.275), (t, Inf, ceiling, 0, 0, 0))
# The ePRC of G.811.1 and the ePRTC of G.8272.1 share their MTIE up to 100 s.
_eprc_mtie(tail::Tuple...) = _pieces((0.1, 1, 4, 0, 0, 0), (1, 100, 3.89, 0, 0, 0.11114), tail...)

const _MASK_TABLE = Pair{Symbol,_Mask}[
    :G811 => _Mask("ITU-T G.811 primary reference clock", "G.811 (09/1997) with Amd. 1 (04/2016), clause 6.1",
        _G811_MTIE, _pieces((0.1, 100, 3, 0, 0, 0), (100, 1000, 0, 0, 0, 0.03), (1000, _under(10_000), 30, 0, 0, 0))),
    :G811_1 => _Mask("ITU-T G.811.1 enhanced primary reference clock", "G.811.1 (08/2017), Tables 1 and 2",
        _eprc_mtie((100, 1000, 15, 0, 0, 0.0375e-3), (1000, Inf, 14.0375, 0, 0, 1e-3)), _pieces((0.1, 10_000, 1, 0, 0, 0))),
    :G812_I => _Mask("ITU-T G.812 node clock, type I", "G.812 (06/2004), Tables 3 and 6",
        _pieces((0.1, 9, 24, 0, 0, 0), (9, 400, 0, 8, 0.5, 0), (400, 10_000, 160, 0, 0, 0)),
        _pieces((0.1, 25, 3, 0, 0, 0), (25, 100, 0, 0, 0, 0.12), (100, 10_000, 12, 0, 0, 0))),
    :G812_II_III => _Mask("ITU-T G.812 node clock, types II and III", "G.812 (06/2004), Tables 4 and 7",
        _pieces((0.1, 1, 40, 0, 0, 0), (1, 10, 0, 40, 0.4, 0), (10, Inf, 100, 0, 0, 0)),
        _pieces((0.1, 2.5, 0, 3.2, -0.5, 0), (2.5, 40, 2, 0, 0, 0), (40, 1000, 0, 0.32, 0.5, 0), (1000, Inf, 10, 0, 0, 0))),
    :G813_opt1 => _Mask("ITU-T G.813 SDH equipment clock, option 1", "G.813 (03/2003), Tables 1 and 3", _EC1_MTIE, _EC1_TDEV),
    :G813_opt2 => _Mask("ITU-T G.813 SDH equipment clock, option 2", "G.813 (03/2003), Tables 4 and 5", _EC2_MTIE, _EC2_TDEV),
    :G8262_opt1 => _Mask("ITU-T G.8262 synchronous Ethernet equipment clock, option 1", "G.8262 (10/2024), Tables 1 and 3",
        _EC1_MTIE, _EC1_TDEV),
    :G8262_opt2 => _Mask("ITU-T G.8262 synchronous Ethernet equipment clock, option 2", "G.8262 (10/2024), Tables 4 and 5",
        _EC2_MTIE, _EC2_TDEV),
    :G8262_1 => _Mask("ITU-T G.8262.1 enhanced synchronous equipment clock", "G.8262.1 (11/2025), Tables 1 and 3",
        _pieces((_under(0.1), 1, 0, 10, 0.155, 0), (1, 100, 0, 10, 0.1, 0), (100, 1000, 0, 6.3, 0.2, 0)),
        _pieces((_under(0.1), 25, 0.64, 0, 0, 0), (25, 100, 0, 0.128, 0.5, 0), (100, 1000, 1.28, 0, 0, 0))),
    :G8272_PRTC_A => _Mask("ITU-T G.8272 primary reference time clock, PRTC-A", "G.8272 (07/2025), Tables 1 and 3",
        _prtc_mtie(273, 100),
        _pieces((0.1, 100, 3, 0, 0, 0), (100, 1000, 0, 0, 0, 0.03), (1000, _under(10_000), 30, 0, 0, 0))),
    :G8272_PRTC_B => _Mask("ITU-T G.8272 primary reference time clock, PRTC-B", "G.8272 (07/2025), Tables 2 and 4",
        _prtc_mtie(54.5, 40),
        _pieces((0.1, 100, 1, 0, 0, 0), (100, 500, 0, 0, 0, 0.01), (500, _under(100_000), 5, 0, 0, 0))),
    :G8272_1 => _Mask("ITU-T G.8272.1 enhanced primary reference time clock",
        "G.8272.1 (01/2024) with Amd. 1 (07/2025), Tables 1 and 2",
        _eprc_mtie((100, 400_000, 15, 0, 0, 0.0375e-3), (400_000, Inf, 30, 0, 0, 0)),
        _pieces((0.1, 30_000, 1, 0, 0, 0), (30_000, 300_000, 0, 0, 0, 3.33333e-5), (300_000, _under(1e6), 10, 0, 0, 0))),
    :G823_PRC => _Mask("ITU-T G.823 network limit at PRC interfaces", "G.823 (03/2000), Tables 6 and 7", _G811_MTIE,
        _pieces((0.1, 100, 3, 0, 0, 0), (100, 1000, 0, 0, 0, 0.03), (1000, 10_000, 30, 0, 0, 0), (10_000, 1e6, 27, 0, 0, 0.0003))),
    :G823_SSU => _Mask("ITU-T G.823 network limit at SSU interfaces", "G.823 (03/2000), Tables 8 and 9",
        _g823_mtie(2.5, 25, 10, 200), _g823_tdev(4.3, 3)),
    :G823_SEC => _Mask("ITU-T G.823 network limit at SEC interfaces", "G.823 (03/2000), Tables 10 and 11",
        _g823_mtie(2.5, 250, 100, 20), _g823_tdev(17.14, 12)),
    :G823_PDH => _Mask("ITU-T G.823 network limit at PDH synchronization interfaces", "G.823 (03/2000), Tables 12 and 13",
        _g823_mtie(7.3, 732, 100, 20), _g823_tdev(48, 34)),
]
const _MASKS = Dict(_MASK_TABLE)

"""
    TELECOM_MASKS

Names of the ITU-T wander masks that [`telecom_mask`](@ref) and
[`mask_check`](@ref) know. Each limits MTIE and TDEV. The first twelve
limit the wander that a clock generates in locked mode, the last four the
wander of a network at an interface.

| Name | Clock or interface | Source |
|:-----|:-------------------|:-------|
| `:G811` | Primary reference clock | [G.811, clause 6.1](@cite ITU-G811) |
| `:G811_1` | Enhanced primary reference clock | [G.811.1, Tables 1 and 2](@cite ITU-G811-1) |
| `:G812_I`, `:G812_II_III` | Node clock of type I, and of types II and III, at constant temperature | [G.812, Tables 3, 4, 6 and 7](@cite ITU-G812) |
| `:G813_opt1`, `:G813_opt2` | SDH equipment clock at constant temperature | [G.813, Tables 1 and 3 to 5](@cite ITU-G813) |
| `:G8262_opt1`, `:G8262_opt2` | Synchronous equipment clock at constant temperature, with the limits of G.813 | [G.8262, Tables 1 and 3 to 5](@cite ITU-G8262) |
| `:G8262_1` | Enhanced synchronous equipment clock at constant temperature | [G.8262.1, Tables 1 and 3](@cite ITU-G8262-1) |
| `:G8272_PRTC_A`, `:G8272_PRTC_B` | Primary reference time clock of class A and of class B | [G.8272, Tables 1 to 4](@cite ITU-G8272) |
| `:G8272_1` | Enhanced primary reference time clock | [G.8272.1, Tables 1 and 2](@cite ITU-G8272-1) |
| `:G823_PRC`, `:G823_SSU`, `:G823_SEC`, `:G823_PDH` | PRC, SSU, SEC and PDH synchronization interfaces | [G.823, Tables 6 to 13](@cite ITU-G823) |

Option 1 of G.813 and G.8262 is for networks of the 2048 kbit/s
hierarchy, option 2 for the 1544 kbit/s hierarchy.

!!! note "The measurement that the recommendations prescribe"
    The recommendations measure MTIE and TDEV through a 10 Hz first-order
    low-pass filter at a sample interval of 1/30 s or less, and TDEV over
    a record of at least ``12τ``. G.8272 and G.8272.1 take a 1 PPS output
    at one sample per second without a filter, where their masks hold
    from ``τ = 1`` s, and a PTP output through a moving average of at
    least 100 time error samples. This package applies no such filter:
    the record is compared as it was measured.

!!! note "Editions"
    The limits are those of the editions in the reference list. G.812 is
    the edition of 2004: its revision of August 2026, and the amendments
    of that month to G.811.1, G.8262.1, G.8272 and G.8272.1, were approved
    but not yet published when the limits were read.
"""
const TELECOM_MASKS = Tuple(first.(_MASK_TABLE))

function _mask_pieces(mask::Symbol, stat::Symbol)
    m = get(_MASKS, mask, nothing)
    m === nothing && throw(ArgumentError("unknown telecom mask :" * String(mask) *
                                         "; use one of " * join(map(String, TELECOM_MASKS), ", ")))
    stat === :mtie && return m.mtie
    stat === :tdev && return m.tdev
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
