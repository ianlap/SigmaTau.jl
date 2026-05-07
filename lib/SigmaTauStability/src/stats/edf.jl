# stats/edf.jl — Equivalent Degrees of Freedom and Confidence Intervals
using Statistics
using Distributions

"""
    calculate_edf(method::Symbol, devs::Vector{Float64}, noises::Vector{Symbol}, m_values::Vector{Int}, taus::Vector{Float64}, N::Int, T::Float64) → Vector{Float64}

Calculates EDF based on Greenhall/Riley approximation tables.
"""
function calculate_edf(method::Symbol, devs::Vector{Float64}, noises::Vector{Symbol}, m_values::Vector{Int}, taus::Vector{Float64}, N::Int, T::Float64)
    edfs = Vector{Float64}(undef, length(devs))
    
    for k in 1:length(devs)
        m = m_values[k]
        noise = noises[k]
        alpha = _alpha_from_noise(noise)
        tau = taus[k]
        
        if method == :adev
            edfs[k] = _calc_edf_core(alpha, 2, m, m, 1, N)
        elseif method == :mdev
            edfs[k] = _calc_edf_core(alpha, 2, m, 1, 1, N)
        elseif method == :hdev
            edfs[k] = _calc_edf_core(alpha, 3, m, m, 1, N)
        elseif method == :mhdev
            edfs[k] = _calc_edf_core(alpha, 3, m, 1, 1, N)
        elseif method == :totdev
            b, c = _coeff_totvar(alpha)
            edfs[k] = b * (T / tau) - c
        elseif method == :mtotdev
            b, c = _coeff_mtot(alpha)
            edfs[k] = b * (T / tau) - c
        elseif method == :htotdev
            b0, b1 = _coeff_htot(alpha)
            edfs[k] = (T / tau) / (b0 + b1 * (tau / T))
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

    sz0 = Float64(_compute_sz(0.0, F, alpha, d))
    bsum = sz0^2
    for j in 1:(J-1)
        szj = _compute_sz(j/S, F, alpha, d)
        bsum += 2 * (1 - j/M) * szj^2
    end
    if J <= M
        szJ = _compute_sz(J/S, F, alpha, d)
        bsum += (1 - J/M) * szJ^2
    end

    bsum > 0 || return NaN
    return M * sz0^2 / bsum
end

function _compute_sw(t::Real, alpha::Int)
    ta = abs(t)
    if alpha == 2;  return -ta
    elseif alpha == 1;  return t^2 * log(max(ta, eps()))
    elseif alpha == 0;  return ta^3
    elseif alpha == -1; return -t^4 * log(max(ta, eps()))
    elseif alpha == -2; return -ta^5
    elseif alpha == -3; return  t^6 * log(max(ta, eps()))
    elseif alpha == -4; return  ta^7
    else; return NaN
    end
end

function _compute_sx(t::Real, F::Int, alpha::Int)
    if F > 100 && alpha <= 0
        return _compute_sw(t, alpha + 2)
    end
    return F^2 * (2 * _compute_sw(t, alpha) -
                  _compute_sw(t - 1/F, alpha) -
                  _compute_sw(t + 1/F, alpha))
end

function _compute_sz(t::Real, F::Int, alpha::Int, d::Int)
    sx(u) = _compute_sx(u, F, alpha)
    if d == 1
        return 2sx(t) - sx(t-1) - sx(t+1)
    elseif d == 2
        return 6sx(t) - 4sx(t-1) - 4sx(t+1) + sx(t-2) + sx(t+2)
    elseif d == 3
        return 20sx(t) - 15sx(t-1) - 15sx(t+1) + 6sx(t-2) + 6sx(t+2) - sx(t-3) - sx(t+3)
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
    alpha == 2  && return (1.90, 2.10)
    alpha == 1  && return (1.20, 1.40)
    alpha == 0  && return (1.10, 1.20)
    alpha == -1 && return (0.85, 0.50)
    alpha == -2 && return (0.75, 0.31)
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
    alpha == 0  && return (0.546, 1.41)
    alpha == -1 && return (0.667, 2.00)
    alpha == -2 && return (0.909, 1.00)
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

Bias factor B(α). Divide raw deviation by B to get unbiased estimate.
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
            table = Dict(0=>-0.005, -1=>-0.149, -2=>-0.229, -3=>-0.283, -4=>-0.321)
            B[k] = 1 / (1 + get(table, clamp(alpha, -4, 0), 0.0))
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
