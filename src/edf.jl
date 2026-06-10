# stats/edf.jl — Equivalent Degrees of Freedom and Confidence Intervals
using Statistics
using Distributions

"""
    calculate_edf(method::Symbol, devs::Vector{Float64}, noises::Vector{Symbol}, m_values::Vector{Int}, taus::Vector{Float64}, N::Int, T::Float64) → Vector{Float64}

Equivalent degrees of freedom for the chosen `method`. Uses the
Greenhall–Riley spectral integrals for ADEV / MDEV / HDEV / MHDEV (via
`_calc_edf_core`), the published EDF coefficient tables for the
total-family estimators, and the SP1065 Table 3 empirical formulas for
`:theo1` (via `_theo1_edf`).

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
            # The SP1065 Table 8 coefficient formula is valid only for
            # τ ≥ 16τ0. Below that floor MTOT reduces to MDEV, so it inherits
            # the modified-Allan EDF (Greenhall–Riley, d=2, F=1) — verified to
            # match Stable32's reported Chi Square DF to ~5 sig figs at small m
            # (and avoids the b·(T/τ) form returning an unphysical EDF > N).
            if m < 16
                edfs[k] = _calc_edf_core(alpha, 2, m, 1, m, N)
            else
                b, c = _coeff_mtot(alpha)
                edfs[k] = b * (T / tau) - c
            end
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
        elseif method == :theo1
            edfs[k] = _theo1_edf(alpha, m, N)
        elseif method == :pdev
            edfs[k] = _pvar_edf(alpha, m, N)
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
    # edf = b·(T/τ) − c, coefficients per noise type from NIST SP1065
    # (Riley 2008) §5.4.3 Table 8 — the published standard, traceable to
    # Greenhall's degrees-of-freedom recipes (Greenhall 1991 / Greenhall &
    # Riley 2003). Used verbatim for all five noise types, matching how
    # `_coeff_totvar` follows SP1065 Table 7.
    alpha == 2  && return (1.90, 2.10)
    alpha == 1  && return (1.20, 1.40)
    alpha == 0  && return (1.10, 1.20)
    alpha == -1 && return (0.85, 0.50)
    alpha == -2 && return (0.75, 0.31)
    return (NaN, NaN)
end

function _coeff_mhtot(alpha::Int)
    # edf = b·(T/τ) − c, valid for τ/τ0 ≥ 16. Measured by Monte Carlo
    # (tools/mc_mhtotdev.jl, full sweep at git e9d209d, seed 20260531, R=3000,
    # N up to 32768) against the shipped mhtotdev pipeline INCLUDING its
    # default global least-squares drift removal (remove_drift=true), so the
    # calibration matches what users compute; all fits R² ≥ 0.9983. MHTOTDEV
    # is novel to SigmaTau, so these have no external reference — see docs
    # theory/mhtotdev_bias_edf.md. Both the Mod-Totvar form below and the
    # TotHvar form (T/τ)/(b0+b1·τ/T) were fit; they tie at the α = ±2 ends
    # (|ΔR²| ≤ 7e-5) and Mod-Totvar fits better elsewhere (up to ΔR² = 2e-3),
    # so the Mod-Totvar form is used uniformly.
    alpha == 2  && return (1.8432, 5.5387)
    alpha == 1  && return (1.2263, 3.9300)
    alpha == 0  && return (1.1129, 3.9205)
    alpha == -1 && return (1.0312, 3.9211)
    alpha == -2 && return (0.8034, 3.0889)
    return (NaN, NaN)
end

# A(α) coefficient of the PVAR EDF approximation, Vernotte–Chen–Rubiola 2020
# (arXiv:2005.13631) Eq. (18): the 3rd-order polynomial
#   A(α) = 27 + (1/4)α + (5/14)α² − (3/4)α³
# fit to massive Monte-Carlo, with A(+2)=23 fixed by the white-PM closed form of
# Vernotte et al. 2016 (arXiv:1506.00687) Eq. (24). Stored as the integer-α values
# the paper reports (the cubic rounded), since noise identification only yields
# integer α. The companion coefficient B(α) = 12 is constant for all α (Eq. 18).
function _pvar_A(alpha::Int)
    alpha == 2  && return 23.0   # cubic 22.93; fixed exact by 2016 Eq. (24)
    alpha == 1  && return 27.0   # cubic 26.86
    alpha == 0  && return 27.0
    alpha == -1 && return 28.0   # cubic 27.86
    alpha == -2 && return 34.0   # cubic 33.93
    return NaN
end

"""
    _pvar_edf(alpha::Int, m::Int, N::Int) → Float64

Equivalent degrees of freedom of PVAR (parabolic variance) for a record of `N`
phase samples at averaging factor `m`, per Vernotte–Chen–Rubiola 2020
(arXiv:2005.13631). The EDF feeds the χ² confidence interval the same way as the
overlapped families (`ν = 2·E[σ²]²/V[σ²]`, their Eq. 19 ≡ the Wave-4 estimator).

Three regimes:

- `m = 1`: PVAR(τ₀) ≡ AVAR(τ₀) (`_pdev_core` returns ADEV at m=1), so the EDF is
  ADEV's Greenhall–Riley value `_calc_edf_core(α, 2, 1, 1, 1, N)`. This also
  sidesteps the paper's note that the closed-form approximation is poor at m=1,2.
- `1 < m ≤ m₁`: the Eq. (16)+(18) approximation
  `ν = 35 / (A(α)·(m/M) − 12·(m/M)²)` with `M = N − 2m`, valid for `m ≲ N/4`.
- `m₁ < m`: the Eq. (22)–(23) semi-log interpolation `ν = a·ln(m) + b` between
  `m₁ = round(1.11·N/4)` (anchored at the Eq.-16 value) and `m₂ = round(0.901·N/2)`
  (anchored at `ν = 1`, the single-estimate limit at `m = N/2`).

`A(α) = 27 + α/4 + 5α²/14 − 3α³/4` (`_pvar_A`); `B = 12`. PVAR is an unbiased
wavelet variance, so no bias correction is applied.
"""
function _pvar_edf(alpha::Int, m::Int, N::Int)
    m <= 0 && return NaN
    m == 1 && return _calc_edf_core(alpha, 2, 1, 1, 1, N)

    A = _pvar_A(alpha)
    isnan(A) && return NaN

    # Eq. (16): closed form, valid for m ≲ N/4.
    edf16(mm) = begin
        M = N - 2 * mm
        M < 1 && return NaN
        r = mm / M
        denom = A * r - 12.0 * r * r
        denom > 0 ? 35.0 / denom : NaN
    end

    m1 = round(Int, 1.11 * N / 4)
    if m <= m1
        return edf16(m)
    end

    # Eqs. (22)–(23): semi-log interpolation across N/4 < m ≤ N/2, with ν=1 at
    # m₂ ≈ N/2 (one effective estimate). Anchored at the Eq.-16 value at m₁.
    m2 = round(Int, 0.901 * N / 2)
    (m1 < 2 || m2 <= m1) && return edf16(m)   # degenerate for tiny N: fall back
    nu_m1 = edf16(m1)
    isfinite(nu_m1) || return NaN
    a = (nu_m1 - 1.0) / (log(m1) - log(m2))
    b = 1.0 - a * log(m2)
    # ν falls to 1 at m₂ ≈ N/2 (single estimate) and is held there for any larger
    # m, rather than extrapolating the line below 1 (Eq. 22–23 / Fig. 2).
    return max(a * log(m) + b, 1.0)
end

"""
    _theo1_edf(alpha::Int, m::Int, N::Int) → Float64

Equivalent χ² degrees of freedom of the Thêo1 statistic at averaging factor
`m` (even) on `N` phase samples, per the empirical formulas of NIST SP1065
(Riley 2008) §5.2.15 Table 3 (Howe & Peppler 2003), with `r = 0.75·m` and the
stated validity condition `τ0 ≤ T/10`. One closed form per power-law noise
type; α outside {−2…2} returns NaN. Also used for ThêoBR / the Thêo1 segment
of ThêoH — the bias-removal ratio is a scalar rescaling, which standard
practice (Stable32) treats as not changing the degrees of freedom.
"""
function _theo1_edf(alpha::Int, m::Int, N::Int)
    m < 2 && return NaN
    r = 0.75 * m
    Nf = Float64(N)
    if alpha == 2        # W PM
        return (0.86 * (Nf + 1.0) * (Nf - 4.0 * r / 3.0) / (Nf - r)) *
               (r / (r + 1.14))
    elseif alpha == 1    # F PM
        return ((4.798 * Nf^2 - 6.374 * Nf * r + 12.387 * r) /
                (sqrt(r + 36.6) * (Nf - r))) * (r / (r + 0.3))
    elseif alpha == 0    # W FM
        return ((4.1 * Nf + 0.8) / r - (3.1 * Nf + 6.5) / Nf) *
               (r^1.5 / (r^1.5 + 5.2))
    elseif alpha == -1   # F FM
        return ((2.0 * Nf^2 - 1.3 * Nf * r - 3.5 * r) / (Nf * r)) *
               (r^3 / (r^3 + 2.3))
    elseif alpha == -2   # RW FM
        return ((4.4 * Nf - 2.0) / (2.9 * r)) *
               (((4.4 * Nf - 1.0)^2 - 8.6 * r * (4.4 * Nf - 1.0) + 11.4 * r^2) /
                (4.4 * Nf - 3.0)^2)
    end
    return NaN
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
- `:mhtot`  — measured by Monte Carlo (`tools/mc_mhtotdev.jl`) over the
  `τ/τ0 ≥ 16` validity window, since MHTOTDEV is novel to SigmaTau and no
  external reference exists. Modeled τ/T-linearly, `B = b0 + b1·(τ/T)` (like
  `:totvar`): `b0 ∈ {1.064, 0.984, 1.019, 1.213, 1.943}` and
  `b1 ∈ {0.017, 0.036, -0.048, -0.321, -3.588}` for α ∈ {2, 1, 0, -1, -2}.
  MHTOTDEV is ≈ unbiased for white/flicker noise (b1 ≈ 0, B ≈ 1) and reads high
  for redder FM (RWFM B ≈ 1.9 at small τ, falling toward 1 as τ→T).
  See docs `theory/mhtotdev_bias_edf.md`.

Anything unrecognised falls through to B = 1.
"""
function bias_correction(noises::Vector{Symbol}, var_type::Symbol, taus::Vector{Float64}, T::Float64)
    B = ones(Float64, length(noises))

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
        elseif var_type == :mhtot
            # B = E[MHTOTVAR]/E[MHVAR] measured by Monte Carlo (tools/mc_mhtotdev.jl,
            # full sweep at git e9d209d) over the τ/τ0 ≥ 16 validity window
            # (Howe et al. 2000), with the wrapper's default global drift removal
            # (remove_drift=true) in the loop. Modeled as B = b0 + b1·(τ/T) — the
            # same τ/T-linear form as :totvar. For the PM noises b1 ≈ 0 (B is
            # constant); for the FM noises the b1 term carries both the boundary
            # bias and the drift-removal suppression near τ → T and cuts the fit
            # residual by ~53 % (WHFM) to ~97 % (RWFM). MHTOTDEV is ≈ unbiased
            # for white/flicker PM (B ≈ 1) and reads high for redder FM
            # (RWFM B ≈ 1.9 at small τ). MHTOTDEV is novel to SigmaTau —
            # no external reference; see docs theory/mhtotdev_bias_edf.md.
            b0_t = Dict(2=>1.0645, 1=>0.9839, 0=>1.0197, -1=>1.2142, -2=>1.9499)
            b1_t = Dict(2=>-0.0034, 1=>-0.0048, 0=>-0.1772, -1=>-0.7024, -2=>-4.9096)
            aa = clamp(alpha, -2, 2)
            B[k] = get(b0_t, aa, 1.0) + get(b1_t, aa, 0.0) * (tau / T)
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
            # A deviation is non-negative by construction, so floor the lower
            # limit at zero — the symmetric Gaussian approximation can otherwise
            # push it below zero for short records at high confidence.
            lower[k] = max(d - half, 0.0)
            upper[k] = d + half
        end
    end

    return lower, upper
end
