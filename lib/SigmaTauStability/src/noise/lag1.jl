# noise/lag1.jl — Lag-1 Autocorrelation for noise identification
using Statistics

const NEFF_RELIABLE = 30

"""
    identify_noise(x::Vector{Float64}, m_values::Vector{Int}; dmin::Int=0, dmax::Int=2) → Vector{Symbol}

Identifies dominant power-law noise type using lag-1 autocorrelation for phase data with B1-ratio fallback.
"""
function identify_noise(x::Vector{Float64}, m_values::Vector{Int}; dmin::Int=0, dmax::Int=2)
    x_clean = _preprocess(x)
    N = length(x_clean)
    noises = Vector{Symbol}(undef, length(m_values))
    last_reliable = :unknown

    for (k, m) in enumerate(m_values)
        N_eff = N ÷ m
        alpha = NaN
        
        try
            if N_eff >= NEFF_RELIABLE
                alpha, _, _, _ = _noise_id_lag1acf(x_clean, m, dmin, dmax)
            else
                alpha, _, _ = _noise_id_b1rn(x_clean, m)
            end
        catch
            alpha = NaN
        end
        
        if !isnan(alpha)
            # Map alpha to symbol
            if round(alpha) == 2
                noises[k] = :WHPM
            elseif round(alpha) == 1
                noises[k] = :FLPM
            elseif round(alpha) == 0
                noises[k] = :WHFM
            elseif round(alpha) == -1
                noises[k] = :FLFM
            elseif round(alpha) == -2
                noises[k] = :RWFM
            else
                noises[k] = :unknown
            end
            last_reliable = noises[k]
        else
            noises[k] = last_reliable
        end
    end

    return noises
end

function _preprocess(x::Vector{Float64})
    x_mean = mean(x)
    x_std  = std(x)
    if x_std < eps()
        return _detrend_linear(x)
    end
    z = abs.((x .- x_mean) ./ x_std)
    return _detrend_linear(x[z .< 5.0])
end

function _detrend_linear(x::Vector{Float64})
    N = length(x)
    N_float = Float64(N)
    sum_i = (N_float * (N_float + 1.0)) / 2.0
    sum_i2 = (N_float * (N_float + 1.0) * (2.0*N_float + 1.0)) / 6.0
    delta = N_float * sum_i2 - sum_i^2
    sum_x = sum(x)
    sum_ix = sum(i * x[i] for i in 1:N)
    
    a = (sum_x * sum_i2 - sum_ix * sum_i) / delta
    b = (N_float * sum_ix - sum_x * sum_i) / delta
    
    return [x[i] - (a + b * i) for i in 1:N]
end

function _detrend_quadratic(x::Vector{Float64})
    N = length(x)
    N_float = Float64(N)
    
    X1 = sum(x)
    X2 = sum(i * x[i] for i in 1:N)
    X3 = sum(i^2 * x[i] for i in 1:N)
    
    S1 = N_float
    S2 = (N_float * (N_float + 1.0)) / 2.0
    S3 = (N_float * (N_float + 1.0) * (2.0*N_float + 1.0)) / 6.0
    S4 = (N_float^2 * (N_float + 1.0)^2) / 4.0
    S5 = (N_float * (N_float + 1.0) * (2.0*N_float + 1.0) * (3.0*N_float^2 + 3.0*N_float - 1.0)) / 30.0
    
    M = [S1 S2 S3; S2 S3 S4; S3 S4 S5]
    V = [X1, X2, X3]
    
    C = M \ V
    a, b, c = C[1], C[2], C[3]
    
    return [x[i] - (a + b*i + c*i^2) for i in 1:N]
end

function _noise_id_lag1acf(x::Vector{Float64}, m::Int, dmin::Int = 0, dmax::Int = 2)
    x_dec = m > 1 ? x[1:m:end] : x
    x_det = _detrend_quadratic(x_dec)
    
    d = 0
    while true
        r1  = _lag1_acf(x_det)
        rho = r1 / (1 + r1)

        if d >= dmin && (rho < 0.25 || d >= dmax)
            p       = -2 * (rho + d)
            alpha   = p + 2
            isnan(alpha) && return (NaN, 0, d, rho)
            return (alpha, round(Int, alpha), d, rho)
        end

        x_det = diff(x_det)
        d += 1
        length(x_det) >= 5 || throw(ArgumentError("Data too short after differencing"))
    end
end

function _lag1_acf(x::Vector{Float64})
    xm   = x .- mean(x)
    ssx  = sum(abs2, xm)
    ssx < eps(Float64) * length(x) && return NaN
    return sum(@view(xm[1:end-1]) .* @view(xm[2:end])) / ssx
end

function _noise_id_b1rn(x::Vector{Float64}, m::Int)
    x_dec = x[1:m:end]
    x_dec = _detrend_quadratic(x_dec)

    avar_val = _simple_avar(x_dec, 1) / Float64(m)^2
    N_avar   = length(x_dec) - 2

    dx = diff(x)
    Nd = (length(dx) ÷ m) * m
    Nd < m && return (NaN, -2, NaN)
    
    y_blocks = reshape(dx[1:Nd], m, :)
    y_avg = vec(mean(y_blocks, dims=1))
    var_class = var(y_avg; corrected=false)

    (isnan(avar_val) || avar_val <= 0) && return (NaN, -2, NaN)

    B1_obs = var_class / avar_val

    mu_list    = [1, 0, -1, -2]
    alpha_list = [-2, -1, 0, 2]
    b1_vals    = [_b1_theory(N_avar, mu) for mu in mu_list]

    mu_best   = mu_list[end]
    alpha_int = alpha_list[end]

    for i in 1:(length(mu_list) - 1)
        boundary = sqrt(b1_vals[i] * b1_vals[i+1])
        if B1_obs > boundary
            mu_best   = mu_list[i]
            alpha_int = alpha_list[i]
            break
        end
    end

    if mu_best == -2
        adev_val = sqrt(avar_val)
        mdev_val = _simple_mdev(x, m, 1.0)
        if !isnan(mdev_val) && adev_val > 0
            Rn_obs = (mdev_val / adev_val)^2
            R_hi = _rn_theory(m, 0)
            R_lo = _rn_theory(m, -1)
            alpha_int = Rn_obs > sqrt(R_hi * R_lo) ? 1 : 2
        end
    end

    return (alpha_int, mu_best, B1_obs)
end

function _b1_theory(N::Int, mu::Int)
    if mu == 2;  return N * (N + 1) / 6
    elseif mu == 1;  return N / 2
    elseif mu == 0;  return N * log(N) / (2 * (N - 1) * log(2))
    elseif mu == -1; return 1.0
    elseif mu == -2; return (N^2 - 1) / (1.5 * N * (N - 1))
    else;            return (N * (1 - N^mu)) / (2 * (N - 1) * (1 - 2^mu))
    end
end

function _rn_theory(af::Int, b::Int)
    if b == 0
        return 1.0 / af
    elseif b == -1
        avar = (1.038 + 3 * log(2π * 0.5 * af)) / (4π^2)
        mvar = 3 * log(256 / 27) / (8π^2)
        return mvar / avar
    else
        return 1.0
    end
end

function _simple_avar(x::Vector{Float64}, m::Int)
    N = length(x)
    L = N - 2m
    L <= 0 && return NaN
    d2 = @view(x[1+2m:N]) .- 2 .* @view(x[1+m:N-m]) .+ @view(x[1:L])
    return sum(abs2, d2) / (2.0 * Float64(m)^2 * L)
end

function _simple_mdev(x::Vector{Float64}, m::Int, tau0::Float64)
    N  = length(x)
    Ne = N - 3m + 1
    Ne <= 0 && return NaN
    cs = cumsum(pushfirst!(copy(x), 0.0))
    s1 = @view(cs[1+m:Ne+m])   .- @view(cs[1:Ne])
    s2 = @view(cs[1+2m:Ne+2m]) .- @view(cs[1+m:Ne+m])
    s3 = @view(cs[1+3m:Ne+3m]) .- @view(cs[1+2m:Ne+2m])
    d  = (s3 .- 2 .* s2 .+ s1) ./ m
    return sqrt(sum(abs2, d) / (2.0 * Float64(m)^2 * tau0^2 * Ne))
end
