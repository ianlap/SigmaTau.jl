# io/outliers.jl — Stable32-style MAD outlier detection and removal.
#
# The check function in Stable32 (NIST SP1065 §10.10, "Check" data function)
# flags frequency samples whose absolute deviation from the median exceeds
# nsigma times the median absolute deviation scaled to a Gaussian sigma
# (MAD / 0.6745). Removal replaces the flagged samples with imputed values
# from the Howe gap-filling machinery so the noise character is preserved.

"""
    find_outliers(data::Union{PhaseData,FrequencyData}; nsigma::Real=5) → Vector{Int}

Return the (sorted, unique) indices of outlying samples using the
Stable32-style median-absolute-deviation test [riley-2008-sp1065](@cite):
a sample is an outlier when

```
|yᵢ − median(y)| / (MAD / 0.6745) > nsigma
```

where `MAD = median(|y − median(y)|)` and 0.6745 rescales the MAD into an
estimate of the standard deviation for normally distributed data. `nsigma=5`
is Stable32's default deglitching threshold.

For `FrequencyData` the test runs directly on the samples `y`. For
`PhaseData` it runs on the first differences `Δᵢ = x[i+1] − x[i]` (a phase
glitch appears as a step in frequency), and each flagged difference `Δᵢ`
flags **phase sample `i+1`** — the sample where the step lands. An isolated
phase glitch at `x[j]` therefore produces two flagged differences (`Δⱼ₋₁` up,
`Δⱼ` down) and flags samples `j` and `j+1`.

`NaN` samples (existing gaps) are excluded from the median/MAD estimates and
are never flagged. Returns an empty vector when nothing exceeds the
threshold.

See also [`remove_outliers`](@ref), [`fillgaps`](@ref).
"""
function find_outliers(data::FrequencyData; nsigma::Real=5)
    nsigma > 0 || throw(ArgumentError("find_outliers: nsigma must be positive, got $nsigma"))
    return _mad_outliers(data.y, nsigma)
end

function find_outliers(data::PhaseData; nsigma::Real=5)
    nsigma > 0 || throw(ArgumentError("find_outliers: nsigma must be positive, got $nsigma"))
    # A phase outlier is a step in the first differences; flagged difference
    # Δᵢ = x[i+1] − x[i] implicates the sample the step lands on, i+1.
    flagged_diffs = _mad_outliers(diff(data.x), nsigma)
    return flagged_diffs .+ 1
end

# MAD test core: indices i with |v[i] − median| > nsigma · MAD/0.6745.
# Statistics computed over finite samples only; NaN entries never flag
# (NaN comparisons are false).
function _mad_outliers(v::AbstractVector{<:Real}, nsigma::Real)
    finite = [x for x in v if !isnan(x)]
    isempty(finite) && return Int[]
    med       = median(finite)
    mad       = median(abs.(finite .- med))
    threshold = nsigma * mad / 0.6745
    return [i for i in eachindex(v) if abs(v[i] - med) > threshold]
end

"""
    remove_outliers(data::Union{PhaseData,FrequencyData}; nsigma::Real=5)
        → PhaseData | FrequencyData

Detect outliers with [`find_outliers`](@ref) and return a new record with the
flagged samples replaced by imputed values from [`fillgaps`](@ref) (Howe &
Schlossberger's reflect-and-filter algorithm, which preserves the local noise
character). The `tau0` and `source` fields carry over; the input is left
untouched. When nothing is flagged the input is returned as-is.

This is the Stable32 deglitching idiom [riley-2008-sp1065](@cite): flag
samples more than `nsigma` scaled MADs from the median, gap them, and fill.

# Examples

```julia
fd = read_frequency("counter_log.csv")
fd_clean = remove_outliers(fd)               # 5-sigma default
fd_strict = remove_outliers(fd; nsigma=3)
```

See also [`find_outliers`](@ref), [`fillgaps`](@ref).
"""
function remove_outliers(data::PhaseData; nsigma::Real=5)
    idx = find_outliers(data; nsigma=nsigma)
    isempty(idx) && return data
    x = Vector{Float64}(data.x)
    x[idx] .= NaN
    return fillgaps(PhaseData(x, data.tau0; source=data.source))
end

function remove_outliers(data::FrequencyData; nsigma::Real=5)
    idx = find_outliers(data; nsigma=nsigma)
    isempty(idx) && return data
    y = Vector{Float64}(data.y)
    y[idx] .= NaN
    return fillgaps(FrequencyData(y, data.tau0; source=data.source))
end
