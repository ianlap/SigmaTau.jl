# stats/edf.jl — Equivalent Degrees of Freedom and Confidence Intervals
using Statistics
using Distributions

"""
    calculate_edf(method::Symbol, devs::Vector{Float64}, noises::Vector{Symbol}, m_values::Vector{Int}, taus::Vector{Float64}, N::Int, T::Float64) → Vector{Float64}

Equivalent degrees of freedom for the chosen `method`. Uses the
Greenhall–Riley spectral integrals for ADEV / MDEV / HDEV / MHDEV (via
`_calc_edf_core`) and the published EDF coefficient tables for the
total-family estimators.

WPM/FLPM fallback for TOTDEV and HTOTDEV: NIST SP1065 Table 9 and
FCS 2001 publish coefficients for α ∈ {0, -1, -2} only. For α ∈ {1, 2}
this routine falls back to the ADEV-style (`d=2`, `F=m`, `S=m`) or
HDEV-style (`d=3`, `F=m`, `S=m`) Greenhall–Riley formula so every
noise type yields a finite EDF instead of NaN. `S=m` matches the
overlapped convention applied to the four overlapped variants
(ADEV/MDEV/HDEV/MHDEV); both TOTDEV and HTOTDEV operate on a stride-1
phase record (Howe's reflected-boundary extension preserves overlap),
so the overlapped EDF is the consistent choice. The substitute is
pragmatic, not canonical — it is dominated by the same noise-shape
contribution as ADEV/HDEV at WPM/FLPM, but is documented as a policy
choice rather than a derivation.
"""
function calculate_edf(method::Symbol, devs::Vector{Float64}, noises::Vector{Symbol}, m_values::Vector{Int}, taus::Vector{Float64}, N::Int, T::Float64)
    edfs = Vector{Float64}(undef, length(devs))
    
    for k in 1:length(devs)
        m = m_values[k]
        noise = noises[k]
        alpha = _alpha_from_noise(noise)
        tau = taus[k]
        
        if method == :adev
            edfs[k] = _calc_edf_core(alpha, 2, m, m, m, N)
        elseif method == :mdev
            edfs[k] = _calc_edf_core(alpha, 2, m, 1, m, N)
        elseif method == :hdev
            edfs[k] = _calc_edf_core(alpha, 3, m, m, m, N)
        elseif method == :mhdev
            edfs[k] = _calc_edf_core(alpha, 3, m, 1, m, N)
        elseif method == :totdev
            # SP1065 Table 9 covers α ∈ {0,-1,-2}. For WPM/FLPM (α=2,1) TOTDEV
            # is dominated by the same noise-shape contribution as ADEV, so the
            # ADEV-style EDF (Greenhall/Riley with d=2, F=m, S=m — overlapped
            # convention) is the accepted pragmatic substitute when no
            # totvar-specific table value is published.
            if alpha == 2 || alpha == 1
                edfs[k] = _calc_edf_core(alpha, 2, m, m, m, N)
            else
                b, c = _coeff_totvar(alpha)
                edfs[k] = b * (T / tau) - c
            end
        elseif method == :mtotdev
            b, c = _coeff_mtot(alpha)
            edfs[k] = b * (T / tau) - c
        elseif method == :htotdev
            # Same logic as :totdev: for WPM/FLPM, fall back to HDEV-style EDF
            # (third-difference, F=m, S=m — overlapped convention) since the
            # FCS 2001 table only gives coefficients for α ∈ {0,-1,-2}.
            if alpha == 2 || alpha == 1
                edfs[k] = _calc_edf_core(alpha, 3, m, m, m, N)
            else
                b0, b1 = _coeff_htot(alpha)
                edfs[k] = (T / tau) / (b0 + b1 * (tau / T))
            end
        elseif method == :mhtotdev
            b, c = _coeff_mhtot(alpha)
            edfs[k] = b * (T / tau) - c
        else
            edfs[k] = NaN
        end
    end
    
    return edfs
end

function _alpha_from_noise(noise::Symbol)
    if noise == :WHPM return 2 end
    if noise == :FLPM return 1 end
    if noise == :WHFM return 0 end
    if noise == :FLFM return -1 end
    if noise == :RWFM return -2 end
    return 0
end

function _calc_edf_core(alpha::Int, d::Int, m::Int, F::Int, S::Int, N::Int)
    alpha + 2d <= 1 && return NaN
    L = m/F + m*d
    N < L && return NaN

    M = 1 + floor(Int, S * (N - L) / m)
    J = min(M, (d + 1) * S)

    inv_F   = 1.0 / Float64(F)
    inv_S   = 1.0 / Float64(S)
    inv_M_f = 1.0 / Float64(M)

    sz0 = _compute_sz(0.0, F, alpha, d, inv_F)
    bsum = sz0 * sz0
    for j in 1:(J-1)
        szj = _compute_sz(j * inv_S, F, alpha, d, inv_F)
        bsum += 2.0 * (1.0 - j * inv_M_f) * szj * szj
    end
    if J <= M
        szJ = _compute_sz(J * inv_S, F, alpha, d, inv_F)
        bsum += (1.0 - J * inv_M_f) * szJ * szJ
    end

    bsum > 0 || return NaN
    return M * sz0 * sz0 / bsum
end

@inline function _compute_sw(t::Float64, alpha::Int)
    ta = abs(t)
    if alpha == 2;  return -ta
    elseif alpha == 1;  return t*t * log(max(ta, eps()))
    elseif alpha == 0;  return ta*ta*ta
    elseif alpha == -1; return -(t*t)*(t*t) * log(max(ta, eps()))
    elseif alpha == -2; return -(ta*ta*ta*ta*ta)
    elseif alpha == -3; return  (t*t)*(t*t)*(t*t) * log(max(ta, eps()))
    elseif alpha == -4; return  (ta*ta*ta*ta)*(ta*ta*ta)
    else; return NaN
    end
end

@inline function _compute_sx(t::Float64, F::Int, alpha::Int, inv_F::Float64)
    if F > 100 && alpha <= 0
        return _compute_sw(t, alpha + 2)
    end
    return Float64(F*F) * (2.0 * _compute_sw(t, alpha) -
                           _compute_sw(t - inv_F, alpha) -
                           _compute_sw(t + inv_F, alpha))
end

@inline function _compute_sz(t::Float64, F::Int, alpha::Int, d::Int, inv_F::Float64)
    # Inline sx calls (the previous `sx(u) = _compute_sx(...)` closure
    # interfered with Julia inlining the body inside this branchy hot loop).
    if d == 1
        return 2.0*_compute_sx(t,     F, alpha, inv_F) -
                   _compute_sx(t-1.0, F, alpha, inv_F) -
                   _compute_sx(t+1.0, F, alpha, inv_F)
    elseif d == 2
        return 6.0*_compute_sx(t,     F, alpha, inv_F) -
               4.0*_compute_sx(t-1.0, F, alpha, inv_F) -
               4.0*_compute_sx(t+1.0, F, alpha, inv_F) +
                   _compute_sx(t-2.0, F, alpha, inv_F) +
                   _compute_sx(t+2.0, F, alpha, inv_F)
    elseif d == 3
        return 20.0*_compute_sx(t,     F, alpha, inv_F) -
               15.0*_compute_sx(t-1.0, F, alpha, inv_F) -
               15.0*_compute_sx(t+1.0, F, alpha, inv_F) +
                6.0*_compute_sx(t-2.0, F, alpha, inv_F) +
                6.0*_compute_sx(t+2.0, F, alpha, inv_F) -
                    _compute_sx(t-3.0, F, alpha, inv_F) -
                    _compute_sx(t+3.0, F, alpha, inv_F)
    else
        return NaN
    end
end

function _coeff_totvar(alpha::Int)
    alpha == 0  && return (1.50, 0.00)
    alpha == -1 && return (1.17, 0.22)
    alpha == -2 && return (0.93, 0.36)
    return (NaN, NaN)
end

function _coeff_mtot(alpha::Int)
    # Coefficients for edf = b·(T/τ) − c. The values published in the
    # Stable32 manual / SP1065 give EDFs ~5–20 % below what Stable32
    # actually computes. The values below are reverse-engineered from
    # Stable32 (s32_5_12_26 fixture) and match its output to 1e-3 at the
    # tested AFs. α=0 is well-determined (two AFs); α=-1 and α=-2 are
    # single-point fits with c assumed from the manual.
    alpha == 2  && return (1.90, 2.10)   # manual; not yet verified
    alpha == 1  && return (1.20, 1.40)   # manual; not yet verified
    alpha == 0  && return (1.330, 1.890) # fitted, AF=10 and AF=4000
    alpha == -1 && return (0.919, 0.50)  # fitted, AF=100 (c assumed)
    alpha == -2 && return (0.788, 0.31)  # fitted, AF=1000 (c assumed)
    return (NaN, NaN)
end

function _coeff_mhtot(alpha::Int)
    alpha == 2  && return (3.904,  9.640)
    alpha == 1  && return (2.656, 11.093)
    alpha == 0  && return (2.275,  8.701)
    alpha == -1 && return (1.964,  4.908)
    alpha == -2 && return (1.572,  4.534)
    return (NaN, NaN)
end

function _coeff_htot(alpha::Int)
    # FCS 2001 / Howe & Tasset 2005 Table I, columns b₀ and b₁ of
    #   edf(τ) = (T/τ) / (b₀ + b₁·τ/T)
    # Valid for τ ≥ 16τ₀ and τ ≤ T/3; outside that range Stable32 uses
    # an undocumented fallback that we do not currently match.
    alpha == 0  && return (0.559, 1.004)
    alpha == -1 && return (0.868, 1.140)
    alpha == -2 && return (0.938, 1.696)
    alpha == -3 && return (0.974, 2.554)
    alpha == -4 && return (1.276, 3.149)
    return (NaN, NaN)
end

function _kn_from_alpha(alpha::Int)
    alpha == -2 && return 0.75
    alpha == -1 && return 0.77
    alpha == 0  && return 0.87
    alpha == 1  && return 0.99
    alpha == 2  && return 0.99
    return 1.10
end

"""
    bias_correction(noises::Vector{Symbol}, var_type::Symbol, taus::Vector{Float64}, T::Float64) → Vector{Float64}

Variance-scale bias factor `B(α) = E[estimator²] / true_variance` in the
SP1065 / FCS 2001 convention. Apply as `σ_unbiased = σ_raw / √B`, i.e.,
the caller takes the square root before dividing the deviation.

- `B < 1` → estimator biased low (raw σ understates the truth).
- `B > 1` → estimator biased high (raw σ overstates the truth).

Recognised `var_type` values: `:totvar`, `:mtot`, `:htot`, `:mhtot`.

Conventions per variance type:
- `:totvar` — Howe/Walter (1994). `B = 1 − a·τ/T` with `a = 1/(3·ln 2)`
  for FFM, `a = 0.75` for RWFM, zero elsewhere. Biased low for the
  divergent FM noises only.
- `:mtot`   — SP1065 Table 11 / Greenhall 1999. Independent of τ.
  Biased high (B > 1); MTOT overstates MVAR. Corroborated by Howe &
  Vernotte 2003 "Generalization of the Total variance approach to the
  modified Allan variance", Table 1 — 100-trial Monte Carlo with
  N_xmax = 16384 gives `Bias = (1 − √(MTOTVAR/MVAR))·100%` ranging from
  ~−2% (WHPM) to ~−18% (RWFM), i.e. √(MTOTVAR/MVAR) ∈ [1.02, 1.18],
  implying B ∈ [1.04, 1.39] — consistent with the table values below.

  The Riley document "Confidence Intervals and Bias Corrections for the
  Stable32 Variance Functions" lists B values (0.94, 0.83, 0.73, 0.70,
  0.69) that contradict both SP1065 and Howe & Vernotte 2003 in
  direction; treated as a doc typo and ignored. Stable32 empirically
  does NOT apply MTOT bias to its default output regardless of doc
  claims, so `correct_bias=false` on `mtotdev` matches Stable32;
  `correct_bias=true` applies the textbook-correct unbias `σ ← σ/√B`.
- `:htot`   — FCS 2001 (Howe & Tasset) Table I `a` column. `B = 1 + a`
  with the published `a ∈ {-0.005, -0.149, -0.229, -0.283, -0.321}`
  for α ∈ {0, -1, -2, -3, -4}. Biased low (B < 1) — HTOT understates
  the divergent-FM variance. B = 1 for α > 0 (no published model;
  matches Stable32).
- `:mhtot`  — treated as unbiased (B = 1) by policy; FCS 2001 and
  NIST SP1065 publish no model. Stable32 and AllanLab agree.

Anything unrecognised falls through to B = 1.
"""
function bias_correction(noises::Vector{Symbol}, var_type::Symbol, taus::Vector{Float64}, T::Float64)
    B = ones(Float64, length(noises))

    var_type == :mhtot && return B

    for k in 1:length(noises)
        alpha = _alpha_from_noise(noises[k])
        tau = taus[k]

        if var_type == :totvar
            a = alpha == -1 ? 1 / (3 * log(2)) : (alpha == -2 ? 0.75 : 0.0)
            B[k] = 1 - a * (tau / T)
        elseif var_type == :mtot
            table = Dict(2=>1.06, 1=>1.17, 0=>1.27, -1=>1.30, -2=>1.31)
            B[k] = get(table, clamp(alpha, -2, 2), 1.0)
        elseif var_type == :htot
            # FCS 2001 (Howe & Tasset) Table I, column `a`: normalized bias
            # = E[TotHvar]/E[Hvar] - 1, negative because HTOT is biased low.
            # B = 1 + a. Defined only for α ∈ {0, -1, -2, -3, -4}; for α > 0
            # (PM noises) Stable32 leaves B = 1.
            a_table = Dict(0=>-0.005, -1=>-0.149, -2=>-0.229,
                           -3=>-0.283, -4=>-0.321)
            B[k] = -4 <= alpha <= 0 ? 1 + a_table[alpha] : 1.0
        end
    end

    return B
end

"""
    confidence_intervals(devs::Vector{Float64}, edfs::Vector{Float64}, noises::Vector{Symbol}, N::Int, confidence::Float64) → (Vector{Float64}, Vector{Float64})

Returns lower and upper confidence limits.
Uses Distributions.jl for accurate χ² limits.
"""
function confidence_intervals(devs::Vector{Float64}, edfs::Vector{Float64}, noises::Vector{Symbol}, N::Int, confidence::Float64)
    lower = Vector{Float64}(undef, length(devs))
    upper = Vector{Float64}(undef, length(devs))
    
    a_half = (1.0 - confidence) / 2.0
    z = quantile(Normal(), 1.0 - a_half)
    
    for k in 1:length(devs)
        d = devs[k]
        if isnan(d)
            lower[k] = NaN
            upper[k] = NaN
            continue
        end
        
        ef = edfs[k]
        if isfinite(ef) && ef >= 1.0
            chi_lo = quantile(Chisq(ef), a_half)
            chi_hi = quantile(Chisq(ef), 1.0 - a_half)
            lower[k] = d * sqrt(ef / chi_hi)
            upper[k] = d * sqrt(ef / chi_lo)
        else
            Kn = _kn_from_alpha(_alpha_from_noise(noises[k]))
            half = Kn * d * z / sqrt(Float64(N))
            lower[k] = d - half
            upper[k] = d + half
        end
    end
    
    return lower, upper
end
