# Equivalent degrees of freedom, bias factors and χ² confidence intervals.

"""
    calculate_edf(method, devs, noises, m_values, taus, N, T) → Vector{Float64}

Equivalent degrees of freedom (EDF) of a variance estimate at each averaging
factor: the ``ν`` of the χ² distribution behind
[`confidence_intervals`](@ref). The estimators call it for the interval
they report.

# Arguments

- `method`: the estimator, one of the symbols below. Any other gives `NaN`.
- `devs`: the deviations. Only their number is used.
- `noises`: the noise type at each point (`:WHPM`, `:FLPM`, `:WHFM`,
  `:FLFM`, `:RWFM`, `:FWFM` or `:RRFM`). Any other label counts as white FM.
- `m_values`, `taus`: the averaging factor and averaging time of each point.
- `N`, `T`: the number of phase samples and the record duration
  ``(N - 1)τ_0``.

# Models

- `:adev`, `:mdev`, `:hdev`, `:mhdev`: the full version of the algorithm of
  [Greenhall and Riley (2003), §3.6](@cite greenhall-2003-edf-stability),
  for the overlapping estimator. The result is `NaN` where the variance
  does not converge for the noise type (flicker-walk and random-run FM with
  `:adev` and `:mdev`).
- `:totdev`: ``b\\,T/τ - c`` on white, flicker and random-walk FM
  ([Greenhall et al. (1999), Eq. (26)](@cite greenhall-1999-totvar), with
  the coefficients as rounded in
  [SP1065 Table 7](@cite riley-2008-sp1065)). On white and flicker PM, the
  model of [Total estimators on PM noise](@ref). At ``m = 1``, the EDF of
  `:adev`.
- `:mtotdev`: ``(T/τ)/(b_0 + b_1τ/T)`` from ``τ = 16τ_0``
  ([Howe and Peppler (2001), Table 2](@cite howe-2001-total-estimators)),
  and the EDF of `:mdev` below it. SP1065 Table 8 gives a different
  approximation, which is not used.
- `:htotdev`: the same form on white FM to random-run FM, with the
  coefficients of
  [Howe et al. (2005), Eq. (7) and Table I](@cite howe-2005-tothvar). The
  paper states it for ``τ ≥ 16τ_0``; it is applied here from ``m = 2``. On
  white and flicker PM, the model of [Total estimators on PM noise](@ref).
  At ``m = 1``, the EDF of `:hdev`.
- `:mhtotdev`: the model of [MHTOTDEV bias and EDF](@ref), and the EDF of
  `:mhdev` at ``m = 1``.
- `:theo1`: the empirical formulas of
  [SP1065 §5.2.15, Table 3](@cite riley-2008-sp1065), for white PM to
  random-walk FM.
- `:pdev`: the approximation of
  [Vernotte et al. (2021), Eqs. (16), (18) and (22) to (24)](@cite vernotte-2021-pvar-noninteger),
  for white PM to random-walk FM, and the EDF of `:adev` at ``m = 1``.
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
            # SP1065 Table 7 covers α ∈ {0,-1,-2}. At m = 1 no difference
            # reaches the extension, so TOTVAR is AVAR term for term.
            if m == 1
                edfs[k] = _calc_edf_core(alpha, 2, m, m, m, N)
            elseif alpha >= 1
                edfs[k] = _totvar_pm_edf(alpha, m, N)
            else
                b, c = _coeff_totvar(alpha)
                edfs[k] = b * (T / tau) - c
            end
        elseif method == :mtotdev
            # Howe and Peppler (2001), Table 2, holds from τ = 16τ0: the
            # Hadamard total column two steps up in α. Below that MTOT takes
            # MDEV's EDF (Greenhall and Riley, d = 2, F = 1), as Stable32
            # does; it reads 20 to 50 % low there, on the wide side.
            if m < 16
                edfs[k] = _calc_edf_core(alpha, 2, m, 1, m, N)
            else
                b0, b1 = _coeff_htot(alpha - 2)
                edfs[k] = (T / tau) / (b0 + b1 * (tau / T))
            end
        elseif method == :htotdev
            # At m = 1 the HTOTDEV kernel is overlapping HDEV (SP1065,
            # Stable32 and allantools all report HDEV there), so it takes
            # HDEV's EDF too; the table formula would give an EDF above N.
            if m == 1
                edfs[k] = _calc_edf_core(alpha, 3, m, m, m, N)
            elseif alpha >= 1
                edfs[k] = _htot_pm_edf(alpha, m, N)
            else
                b0, b1 = _coeff_htot(alpha)
                edfs[k] = (T / tau) / (b0 + b1 * (tau / T))
            end
        elseif method == :mhtotdev
            # At m = 1 MHTOTVAR is 124/225 of MHVAR term for term (see `mhtotdev`).
            if m == 1
                edfs[k] = _calc_edf_core(alpha, 3, 1, 1, 1, N)
            else
                b, c, d, e = _coeff_mhtot(alpha)
                s = tau / T
                edfs[k] = b / s - c + d * s + e * s^3
                alpha == 2 && m < 16 && (edfs[k] *= _MHTOT_WHPM_EDF_SMALL[m-1])
            end
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
    if noise == :FWFM return -3 end
    if noise == :RRFM return -4 end
    return 0
end

# Greenhall and Riley (2003), "Uncertainty of stability variances based on
# finite differences", §3: the EDF algorithm that Stable32 and allantools
# implement. `F` is the filter factor (1 for the modified variances, m for
# the unmodified ones) and `S` the stride factor (m for the overlapped
# estimators, 1 for the non-overlapped). This is the full version (§3.6).
# The basic sum (Eq. (10)) is evaluated whenever its length J is at most
# J_MAX = 100. Beyond that, with r = M/S ≥ d + 1, Tables 1 to 3 give the
# large-J asymptote; otherwise the substitution m′ = J_MAX/r shortens the
# sum. For unmodified variances with α ≤ 0 the filter factor is taken as ∞
# once m(d + 1) > J_MAX (Case 2). Following every branch reproduces
# Stable32's χ² bounds. The simplified version (§3.5, Eq. (13)), which sums
# the whole series, differs by up to about 4 % in EDF where J > 100, on
# white FM and flicker PM under the unmodified variances.
const _GR_JMAX = 100

# Table 1 (modified variances, F = 1): (a0, a1) by α = 2 … −4 (rows), d = 1 … 3.
const _GR_TABLE1 = (
    ((2/3, 1/3),       (7/9, 1/2),       (22/25, 2/3)),
    ((0.840, 0.345),   (0.997, 0.616),   (1.141, 0.843)),
    ((1.079, 0.368),   (1.033, 0.607),   (1.184, 0.848)),
    ((NaN, NaN),       (1.048, 0.534),   (1.180, 0.816)),
    ((NaN, NaN),       (1.302, 0.535),   (1.175, 0.777)),
    ((NaN, NaN),       (NaN, NaN),       (1.194, 0.703)),
    ((NaN, NaN),       (NaN, NaN),       (1.489, 0.702)),
)
# Table 2 (unmodified variances, F = m): (a0, a1) by α = 2 … −4, d = 1 … 3.
const _GR_TABLE2 = (
    ((3/2, 1/2),       (35/18, 1.0),     (231/100, 3/2)),
    ((78.6, 25.2),     (790.0, 410.0),   (9950.0, 6520.0)),
    ((2/3, 1/6),       (2/3, 1/3),       (7/9, 1/2)),
    ((NaN, NaN),       (0.852, 0.375),   (0.997, 0.617)),
    ((NaN, NaN),       (1.079, 0.368),   (1.033, 0.607)),
    ((NaN, NaN),       (NaN, NaN),       (1.053, 0.553)),
    ((NaN, NaN),       (NaN, NaN),       (1.302, 0.535)),
)
# Table 3 (unmodified, α = 1): (b0, b1) by d = 1 … 3.
const _GR_TABLE3 = ((6.0, 4.0), (15.23, 12.0), (47.8, 40.0))

_gr_table(t, alpha, d) = t[3 - alpha][d]

function _calc_edf_core(alpha::Int, d::Int, m::Int, F::Int, S::Int, N::Int)
    alpha + 2d <= 1 && return NaN
    L = m/F + m*d
    N < L && return NaN

    M = 1 + floor(Int, S * (N - L) / m)
    J = min(M, (d + 1) * S)
    r = M / S
    Jm = _GR_JMAX

    if F == 1                                    # case 1: modified variances
        J <= Jm && return 1 / _gr_inv_edf(alpha, d, J, M, Float64(S), 1.0)
        if r >= d + 1
            a0, a1 = _gr_table(_GR_TABLE1, alpha, d)
            return 1 / ((a0 - a1 / r) / r)
        end
        return 1 / _gr_inv_edf(alpha, d, Jm, Jm, Jm / r, 1.0)
    end

    F == m || return 1 / _gr_inv_edf(alpha, d, J, M, Float64(S), Float64(F))

    if alpha <= 0                                # case 2: unmodified, α ≤ 0
        if J <= Jm
            Fp = m * (d + 1) <= Jm ? Float64(m) : Inf
            return 1 / _gr_inv_edf(alpha, d, J, M, Float64(S), Fp)
        elseif r >= d + 1
            a0, a1 = _gr_table(_GR_TABLE2, alpha, d)
            return 1 / ((a0 - a1 / r) / r)
        end
        return 1 / _gr_inv_edf(alpha, d, Jm, Jm, Jm / r, Inf)
    elseif alpha == 1                            # case 3: unmodified, flicker PM
        J <= Jm && return 1 / _gr_inv_edf(alpha, d, J, M, Float64(S), Float64(m))
        b0, b1 = _GR_TABLE3[d]
        if r >= d + 1
            a0, a1 = _gr_table(_GR_TABLE2, alpha, d)
            return 1 / ((a0 - a1 / r) / ((b0 + b1 * log(m))^2 * r))
        end
        mp = Jm / r
        return 1 / (_gr_basic_sum(alpha, d, Jm, Jm, mp, mp) / ((b0 + b1 * log(m))^2 * Jm))
    end
    # case 4: unmodified, white PM. Closed form once K = ⌈r⌉ > d; the exact sum
    # (to which the closed form is equal) covers the K ≤ d corner the paper
    # leaves open.
    if ceil(r) > d
        a0 = binomial(4d, 2d) / binomial(2d, d)^2
        return 1 / ((a0 - (d / 2) / r) / M)
    end
    return 1 / _gr_inv_edf(alpha, d, J, M, Float64(S), Float64(m))
end

# 1/edf = BasicSum(J, M, S, F) / (s_z(0)² · M): Greenhall and Riley (2003), Eq. (13).
function _gr_inv_edf(alpha::Int, d::Int, J::Int, M::Int, S::Float64, F::Float64)
    sz0 = _compute_sz(0.0, F, alpha, d)
    return _gr_basic_sum(alpha, d, J, M, S, F) / (sz0 * sz0 * M)
end

# BasicSum of Greenhall and Riley (2003), Eq. (10).
function _gr_basic_sum(alpha::Int, d::Int, J::Int, M::Int, S::Float64, F::Float64)
    sz0 = _compute_sz(0.0, F, alpha, d)
    bsum = sz0 * sz0
    inv_M = 1.0 / M
    for j in 1:(J-1)
        szj = _compute_sz(j / S, F, alpha, d)
        bsum += 2.0 * (1.0 - j * inv_M) * szj * szj
    end
    szJ = _compute_sz(J / S, F, alpha, d)
    bsum += (1.0 - J * inv_M) * szJ * szJ
    return bsum
end

@inline function _compute_sw(t::Float64, alpha::Int)
    ta = abs(t)
    if alpha == 2;  return -ta
    elseif alpha == 1;  return t == 0 ? 0.0 : t*t * log(ta)
    elseif alpha == 0;  return ta*ta*ta
    elseif alpha == -1; return t == 0 ? 0.0 : -(t*t)*(t*t) * log(ta)
    elseif alpha == -2; return -(ta*ta*ta*ta*ta)
    elseif alpha == -3; return t == 0 ? 0.0 : (t*t)*(t*t)*(t*t) * log(ta)
    elseif alpha == -4; return  (ta*ta*ta*ta)*(ta*ta*ta)
    else; return NaN
    end
end

# s_x(t) of Greenhall and Riley (2003), Eq. (8); F = Inf is the limit
# F → ∞, s_w(t, α + 2).
@inline function _compute_sx(t::Float64, F::Float64, alpha::Int)
    isinf(F) && return _compute_sw(t, alpha + 2)
    inv_F = 1.0 / F
    return F*F * (2.0 * _compute_sw(t, alpha) -
                  _compute_sw(t - inv_F, alpha) -
                  _compute_sw(t + inv_F, alpha))
end

@inline function _compute_sz(t::Float64, F::Float64, alpha::Int, d::Int)
    if d == 1
        return 2.0*_compute_sx(t,     F, alpha) -
                   _compute_sx(t-1.0, F, alpha) -
                   _compute_sx(t+1.0, F, alpha)
    elseif d == 2
        return 6.0*_compute_sx(t,     F, alpha) -
               4.0*_compute_sx(t-1.0, F, alpha) -
               4.0*_compute_sx(t+1.0, F, alpha) +
                   _compute_sx(t-2.0, F, alpha) +
                   _compute_sx(t+2.0, F, alpha)
    elseif d == 3
        return 20.0*_compute_sx(t,     F, alpha) -
               15.0*_compute_sx(t-1.0, F, alpha) -
               15.0*_compute_sx(t+1.0, F, alpha) +
                6.0*_compute_sx(t-2.0, F, alpha) +
                6.0*_compute_sx(t+2.0, F, alpha) -
                    _compute_sx(t-3.0, F, alpha) -
                    _compute_sx(t+3.0, F, alpha)
    else
        return NaN
    end
end

# b and c of edf = b·T/τ − c for α = 0, −1, −2: SP1065 Table 7, which rounds
# Greenhall et al. (1999), Table I (3/2 and 0; 24(ln 2)²/π² and 0.222;
# 140/151 and 0.358) to two decimals.
function _coeff_totvar(alpha::Int)
    alpha == 0  && return (1.50, 0.00)
    alpha == -1 && return (1.17, 0.22)
    alpha == -2 && return (0.93, 0.36)
    return (NaN, NaN)
end

# edf = b/s − c + d·s + e·s³ with s = τ/T, for α = 2 … −4: fitted to the exact
# EDF of the statistic over 16 ≤ τ/τ0 and τ ≤ T/4 (tools/mhtot_models.jl),
# largest residual 1 to 4.5 %. The estimator is defined in this package, so
# no published table exists; see the "MHTOTDEV bias and EDF" docs page.
function _coeff_mhtot(alpha::Int)
    alpha == 2  && return (1.5949, 6.535, 20.18, -52.1)
    alpha == 1  && return (1.0656, 4.223, 10.19, 7.8)
    alpha == 0  && return (0.9420, 3.506, 6.50, 23.3)
    alpha == -1 && return (0.8938, 3.124, 3.94, 31.2)
    alpha == -2 && return (0.8596, 2.676, 0.43, 48.8)
    alpha == -3 && return (0.8129, 2.115, -4.27, 83.8)
    alpha == -4 && return (0.6508, 1.991, -2.63, 78.3)
    return (NaN, NaN, NaN, NaN)
end

# A(α) of the PVAR EDF approximation: Vernotte, Chen and Rubiola (2021),
# Eq. (18), the cubic A(α) = 27 + α/4 + 5α²/14 − 3α³/4 fitted to Monte Carlo
# simulation, with A(+2) = 23 from the white-PM closed form of Vernotte et
# al. (2016), Eq. (24). The integer-α values the paper reports (the cubic,
# rounded) are stored, since noise identification yields integer α. The
# companion coefficient B = 12 for every α.
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

Equivalent degrees of freedom of the parabolic variance on `N` phase samples
at averaging factor `m`, after Vernotte, Chen and Rubiola (2021):

- `m = 1`: the EDF of the overlapping Allan variance, which the statistic
  equals there. The paper notes that its approximation is poor at m = 1
  and 2.
- `1 < m ≤ m₁`: `ν = 35 / (A·m/M − 12·(m/M)²)` with `M = N − 2m` (their
  Eqs. (16) and (18)), where `m₁ = round(1.11·N/4)`.
- `m > m₁`: `ν = a·ln m + b`, the line in ln m from the value at `m₁` down
  to `ν = 1` at `m₂ = round(0.901·N/2)` (their Eqs. (22) to (24)), and 1
  beyond `m₂`.

`A` is `_pvar_A(alpha)`. Flicker-walk and random-run FM, outside the range
the paper covers, give `NaN`.
"""
function _pvar_edf(alpha::Int, m::Int, N::Int)
    m <= 0 && return NaN
    m == 1 && return _calc_edf_core(alpha, 2, 1, 1, 1, N)

    A = _pvar_A(alpha)
    isnan(A) && return NaN

    # Eq. (16): closed form, for m up to about N/4.
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

    # Eqs. (22) to (24): interpolation in ln m across N/4 < m ≤ N/2, from the
    # Eq. (16) value at m₁ down to ν = 1 at m₂ (one effective estimate).
    m2 = round(Int, 0.901 * N / 2)
    (m1 < 2 || m2 <= m1) && return edf16(m)   # too few samples for the line
    nu_m1 = edf16(m1)
    isfinite(nu_m1) || return NaN
    a = (nu_m1 - 1.0) / (log(m1) - log(m2))
    b = 1.0 - a * log(m2)
    # ν is held at 1 beyond m₂, not extrapolated below it.
    return max(a * log(m) + b, 1.0)
end

"""
    _theo1_edf(alpha::Int, m::Int, N::Int) → Float64

Equivalent degrees of freedom of Thêo1 at the even averaging factor `m` on
`N` phase samples: the empirical formulas of SP1065 §5.2.15, Table 3, with
`r = 0.75·m`, stated there for `τ0 ≤ T/10`. One formula per noise type from
white PM to random-walk FM; any other gives `NaN`. ThêoBR and the Thêo1
segment of ThêoH take the same EDF, as in SP1065 §5.4.4.
"""
function _theo1_edf(alpha::Int, m::Int, N::Int)
    m < 2 && return NaN
    r = 0.75 * m
    Nf = Float64(N)
    if alpha == 2        # WHPM
        return (0.86 * (Nf + 1.0) * (Nf - 4.0 * r / 3.0) / (Nf - r)) *
               (r / (r + 1.14))
    elseif alpha == 1    # FLPM
        return ((4.798 * Nf^2 - 6.374 * Nf * r + 12.387 * r) /
                (sqrt(r + 36.6) * (Nf - r))) * (r / (r + 0.3))
    elseif alpha == 0    # WHFM
        return ((4.1 * Nf + 0.8) / r - (3.1 * Nf + 6.5) / Nf) *
               (r^1.5 / (r^1.5 + 5.2))
    elseif alpha == -1   # FLFM
        return ((2.0 * Nf^2 - 1.3 * Nf * r - 3.5 * r) / (Nf * r)) *
               (r^3 / (r^3 + 2.3))
    elseif alpha == -2   # RWFM
        return ((4.4 * Nf - 2.0) / (2.9 * r)) *
               (((4.4 * Nf - 1.0)^2 - 8.6 * r * (4.4 * Nf - 1.0) + 11.4 * r^2) /
                (4.4 * Nf - 3.0)^2)
    end
    return NaN
end

function _coeff_htot(alpha::Int)
    # b₀ and b₁ of edf(τ) = (T/τ) / (b₀ + b₁·τ/T): Howe and Peppler (2001),
    # Table 3, and Howe et al. (2005), Table I, stated for τ ≥ 16τ₀ and
    # τ ≤ T/3. Stable32 uses a fallback outside that range that its manual
    # does not describe; it is not reproduced here. Table 2 of Howe and
    # Peppler gives the modified total variance the same values at α + 2.
    alpha == 0  && return (0.559, 1.004)
    alpha == -1 && return (0.868, 1.140)
    alpha == -2 && return (0.938, 1.696)
    alpha == -3 && return (0.974, 2.554)
    alpha == -4 && return (1.276, 3.149)
    return (NaN, NaN)
end

# Normalized bias a = E[total]/E[parent] − 1 of the modified total variance
# for α = 2, 1, 0, −1, −2 (Howe and Peppler (2001), Table 2); the Hadamard
# total variance takes the same column at α = 0 … −4 (their Table 3; Howe et
# al. (2005), Table I).
const _MTOT_NBIAS = (-0.005, -0.149, -0.229, -0.283, -0.321)

# The published values hold from τ/τ0 = 16. Below it, the exact expectation
# at m = 2 … 15, same rows (tools/total_pm_models.jl). White PM dips at odd m;
# the flicker rows are the mean of two flicker models, each within 2.4 % of
# it at m = 2 and 1 % from m = 5.
const _MTOT_B_SMALL = (
    (0.9979, 0.9375, 0.9960, 0.9726, 0.9956, 0.9828, 0.9955, 0.9871, 0.9954, 0.9894, 0.9954, 0.9908, 0.9953, 0.9917),
    (0.8772, 0.8435, 0.8612, 0.8520, 0.8566, 0.8525, 0.8546, 0.8524, 0.8536, 0.8522, 0.8529, 0.8520, 0.8525, 0.8518),
    (0.8144, 0.7811, 0.7837, 0.7765, 0.7766, 0.7739, 0.7740, 0.7727, 0.7727, 0.7720, 0.7721, 0.7716, 0.7717, 0.7713),
    (0.7354, 0.7216, 0.7219, 0.7192, 0.7192, 0.7182, 0.7183, 0.7178, 0.7179, 0.7176, 0.7176, 0.7174, 0.7175, 0.7174),
    (0.6866, 0.6814, 0.6812, 0.6800, 0.6799, 0.6794, 0.6794, 0.6791, 0.6792, 0.6790, 0.6790, 0.6789, 0.6790, 0.6789),
)
_mtot_bias(row::Int, m::Int) = m < 16 ? _MTOT_B_SMALL[row][m-1] : 1 + _MTOT_NBIAS[row]

# B = E[MHTOTVAR]/E[MHVAR] for α = 2 … −4: the exact expectation of the
# statistic (tools/mhtot_models.jl). It depends on m only: the large-m value,
# within 0.4 % from τ/τ0 = 16, and the values at m = 2 … 15 (flicker rows the
# mean of two flicker models, each within 1.3 % of it at m = 2).
const _MHTOT_B = (0.9000, 0.8218, 0.7917, 0.7828, 0.7867, 0.7996, 0.8196)
const _MHTOT_B_SMALL = (
    (0.9000, 0.8584, 0.9013, 0.8857, 0.9012, 0.8931, 0.9011, 0.8961, 0.9009, 0.8976, 0.9008, 0.8984, 0.9007, 0.8989),
    (0.8475, 0.8233, 0.8340, 0.8275, 0.8292, 0.8266, 0.8269, 0.8256, 0.8256, 0.8249, 0.8248, 0.8243, 0.8243, 0.8239),
    (0.8281, 0.8059, 0.8064, 0.8013, 0.7999, 0.7980, 0.7972, 0.7962, 0.7957, 0.7952, 0.7948, 0.7945, 0.7943, 0.7940),
    (0.8052, 0.7948, 0.7927, 0.7900, 0.7887, 0.7876, 0.7869, 0.7864, 0.7859, 0.7856, 0.7853, 0.7851, 0.7849, 0.7847),
    (0.8009, 0.7964, 0.7946, 0.7928, 0.7917, 0.7909, 0.7903, 0.7898, 0.7895, 0.7892, 0.7890, 0.7888, 0.7886, 0.7885),
    (0.8025, 0.8063, 0.8057, 0.8046, 0.8039, 0.8033, 0.8028, 0.8024, 0.8021, 0.8019, 0.8017, 0.8015, 0.8014, 0.8013),
    (0.8145, 0.8239, 0.8243, 0.8237, 0.8232, 0.8227, 0.8223, 0.8220, 0.8218, 0.8216, 0.8214, 0.8213, 0.8211, 0.8210),
)

# Exact EDF of MHTOTVAR on white PM over the `_coeff_mhtot` model at
# m = 2 … 15, where the model, fitted from m = 16, reads high.
const _MHTOT_WHPM_EDF_SMALL = (0.654, 0.746, 0.878, 0.889, 0.941, 0.941, 0.968, 0.966, 0.982, 0.980, 0.991, 0.989, 0.997, 0.995)

# The total estimators on white and flicker PM, where nothing is published:
# bias and EDF from the exact moments of each statistic
# (tools/total_pm_models.jl; docs total_pm.md).

# N / EDF of the overlapping Allan variance on flicker PM for N ≫ m, in the
# sampled-process model of Greenhall & Riley: a fit to their algorithm.
_flpm_avar_g(m::Int) = (L = log(m); m / (0.789 + 0.276L + 0.204L^2))

# TOTVAR is AVAR plus the terms that reach the reflection. Each of those
# carries an end sample at weight 2, so the two end samples add a variance
# that grows as m². White PM in closed form: interior covariances 6, −4, 1
# give 140N, the shared end samples 64m², the remaining edge terms 225m.
function _totvar_pm_edf(alpha::Int, m::Int, N::Int)
    B2 = _total_bias(:totvar, alpha, m, N)^2
    alpha == 2 && return 72.0 * N^2 * B2 / (140N + 225m + 64.0 * m^2)
    return B2 / (_flpm_avar_g(m) / N + 0.93 * (m / N)^2)
end

# TOTHVAR on PM: N/g interior terms, thinned as the 3m-sample windows come to
# overlap (levelling at s = τ/T = s0) and again as the last windows run out.
function _htot_pm_edf(alpha::Int, m::Int, N::Int)
    g, c, s0, a, b = alpha == 2 ? (2.45, 2.0, 0.125, 1.0, 0.03) : (_flpm_avar_g(m), 2.5, 0.16, 0.6, 0.04)
    s = m / (N - 1)
    return _total_bias(:htot, alpha, m, N)^2 * (N / g) * (1 - a * exp(-(N - 3m) / (b * N))) / (1 + c * min(s, s0))
end

# Variance-scale bias factor E[total] / E[parent] at averaging factor m on N
# phase samples. `published` keeps to the published tables: nothing on PM
# for TOTVAR and TOTHVAR, and the τ ≥ 16τ0 value of MTOT and TOTHVAR at
# every m > 1.
function _total_bias(var_type::Symbol, alpha::Int, m::Int, N::Int, published::Bool=false)
    if published
        var_type in (:totvar, :htot) && alpha >= 1 && return 1.0
        var_type in (:mtot, :htot) && m > 1 && (m = max(m, 16))
    end
    s = m / (N - 1)
    if var_type == :totvar
        alpha == 2  && return 1 + (8 / 3) * (m ÷ 2) / (N - 2)
        alpha == 1  && return 1 + (0.593 + 0.068 * log(m)) * s
        alpha == -1 && return 1 - s / (3 * log(2))
        alpha == -2 && return 1 - 0.75 * s
    elseif var_type == :mtot
        m == 1 && return 0.5
        alpha >= -2 && return _mtot_bias(3 - alpha, m)
    elseif var_type == :htot
        m == 1 && return 1.0
        alpha == 2 && return isodd(m) ? 1.3697 - 0.8198 / m : 1.3838 - 0.2759 / m
        alpha == 1 && return (5.148 + 1.440 * log(m)) / (4.75 + log(m)) - (isodd(m) ? 0.396 : 0.114) / m
        return _mtot_bias(1 - alpha, m)
    elseif var_type == :mhtot
        return m == 1 ? 124 / 225 : m < 16 ? _MHTOT_B_SMALL[3-alpha][m-1] : _MHTOT_B[3-alpha]
    end
    return 1.0
end

"""
    bias_correction(noises, var_type, m_values, N; published=false) → Vector{Float64}

Bias factor ``B = E[\\text{total variance}]/E[\\text{parent variance}]`` of a
total estimator at each averaging factor, for a record of `N` phase samples.
The estimators divide a total deviation by ``\\sqrt{B}``, so ``B < 1`` means
that the uncorrected deviation reads low. `noises` gives the noise type at
each entry of `m_values`; a label outside the seven power-law names counts
as white FM. Below, ``s = τ/T = m/(N-1)``.

- `var_type = :totvar` (total variance against Allan variance):
  ``B = 1 - a\\,s`` with ``a = 1/(3\\ln 2)`` for flicker FM and ``a = 3/4``
  for random-walk FM
  ([Greenhall et al. (1999), Eq. (25) and Table I](@cite greenhall-1999-totvar);
  [SP1065 §5.11, Eq. (52)](@cite riley-2008-sp1065)), and ``B = 1`` for
  white FM. On PM noise the estimator reads high:
  ``B = 1 + (8/3)⌊m/2⌋/(N-2)`` for white PM, which is exact, and
  ``B = 1 + (0.593 + 0.068\\ln m)\\,s`` for flicker PM.
- `:mtot` (modified total against modified Allan): ``B = 0.995``, 0.851,
  0.771, 0.717 and 0.679 for white PM, flicker PM, white FM, flicker FM and
  random-walk FM, independent of ``τ`` from ``τ = 16τ_0``
  ([Howe and Peppler (2001), Table 2](@cite howe-2001-total-estimators),
  which lists ``a = B - 1``). At ``m = 2, …, 15`` the exact expectation at
  that ``m`` is used, up to 6 % from the published value, and ``B = 1/2``
  at ``m = 1``, exactly.
- `:htot` (Hadamard total against Hadamard): the same five values for white
  FM to random-run FM
  ([Howe et al. (2005), Table I](@cite howe-2005-tothvar)), with the same
  values below ``m = 16``. On PM noise the estimator reads high:
  ``B = 1.384 - 0.28/m`` (even ``m``) or ``1.370 - 0.82/m`` (odd ``m``) for
  white PM, and ``(5.15 + 1.44\\ln m)/(4.75 + \\ln m)`` less ``0.11/m``
  (even) or ``0.40/m`` (odd) for flicker PM. ``B = 1`` at ``m = 1``, where
  the estimate is the Hadamard deviation.
- `:mhtot` (modified Hadamard total against modified Hadamard):
  ``B = 0.900``, 0.822, 0.792, 0.783, 0.787, 0.800 and 0.820 for white PM
  to random-run FM from ``m = 16``, the exact value at each
  ``m = 2, …, 15``, and 124/225 at ``m = 1``.

A noise type not listed for a `var_type`, or any other `var_type`, gives
``B = 1``.

The PM values of `:totvar` and `:htot`, every value below ``m = 16`` and all
of `:mhtot` have no published source: they are the expectation of each
statistic, computed in this package (see
[Total estimators on PM noise](@ref) and [MHTOTDEV bias and EDF](@ref)).
`published=true` leaves them out. It gives ``B = 1`` on PM for `:totvar` and
`:htot`, and the published value of `:mtot` and `:htot` at every ``m > 1``.
`:mhtot` is unchanged.

!!! note "SP1065 Table 11"
    [SP1065 §5.12](@cite riley-2008-sp1065) lists modified total bias
    factors of 1.06 to 1.31 and describes them as the ratio of the expected
    modified total variance to the modified Allan variance. The estimator
    reads low, so that ratio is below 1, and those values are not used.
"""
bias_correction(noises::Vector{Symbol}, var_type::Symbol, m_values::Vector{Int}, N::Int; published::Bool=false) =
    [_total_bias(var_type, _alpha_from_noise(noises[k]), m_values[k], N, published) for k in eachindex(noises)]

"""
    confidence_intervals(devs, edfs, confidence) → (lower, upper)

χ² confidence limits of each deviation ``σ`` from its equivalent degrees of
freedom ``ν`` ([SP1065 §5.3.2](@cite riley-2008-sp1065)):

```math
σ\\sqrt{\\frac{ν}{χ^2_{(1+p)/2}(ν)}} \\quad\\text{and}\\quad σ\\sqrt{\\frac{ν}{χ^2_{(1-p)/2}(ν)}} ,
```

where ``p`` is `confidence` and ``χ^2_q(ν)`` is the ``q`` quantile of the χ²
distribution with ``ν`` degrees of freedom. Every ``ν > 0`` has an
interval, however wide. A limit is `NaN` where the deviation is not finite
or ``ν`` is not positive.
"""
function confidence_intervals(devs::Vector{Float64}, edfs::Vector{Float64}, confidence::Float64)
    a = (1 - confidence) / 2
    lower = fill(NaN, length(devs))
    upper = fill(NaN, length(devs))
    for k in eachindex(devs)
        d, ν = devs[k], edfs[k]
        isfinite(d) && isfinite(ν) && ν > 0 || continue
        lower[k] = d * sqrt(ν / quantile(Chisq(ν), 1 - a))
        upper[k] = d * sqrt(ν / quantile(Chisq(ν), a))
    end
    return lower, upper
end
