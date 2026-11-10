# Outlier detection by the median absolute deviation (SP1065 §10.11; the
# Check function of Stable32): a frequency sample is an outlier when it lies
# more than nsigma scaled MADs from the median, the scaled MAD being
# MAD / 0.6745, the standard deviation of normal data. Removal replaces the
# outliers by filled values (`fillgaps`), where Stable32 leaves gaps.

"""
    find_outliers(data; nsigma=5) → Vector{Int}

Indices of the outlying samples of a record, by the median absolute
deviation test of [SP1065 §10.11](@cite riley-2008-sp1065) and of the Check
function of [Stable32](@cite hamilton-2008-stable32-manual). A fractional
frequency sample ``y_i`` is an outlier when

```math
\\lvert y_i - m \\rvert > n_σ \\frac{\\mathrm{MAD}}{0.6745},
\\qquad m = \\operatorname{median}(y),
\\quad \\mathrm{MAD} = \\operatorname{median} \\lvert y - m \\rvert .
```

The factor 0.6745 makes the scaled MAD equal to the standard deviation of
normally distributed data. The default `nsigma=5` is the limit that SP1065
calls common and the default of Stable32.

For a `FrequencyData` record the test is made on its samples. For a
`PhaseData` record it is made on the first differences, because a phase
jump is a spike in frequency, and a flagged difference ``x_{i+1} - x_i``
marks phase sample ``i + 1``, where the jump lands. A single displaced
phase sample ``x_j`` therefore marks samples ``j`` and ``j + 1``.

Missing samples are left out of the median and of the MAD and are never
marked. The test is applied once; the automatic removal of Stable32 repeats
it until no outlier is left. The result is empty when no sample exceeds the
limit, and when the MAD is zero, as in a coarsely quantized record where
more than half of the samples are equal.
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
    med       = _median_nonan(finite)
    mad       = _median_nonan(abs.(finite .- med))
    mad > 0 || return Int[]
    threshold = nsigma * mad / 0.6745
    return [i for i in eachindex(v) if abs(v[i] - med) > threshold]
end

"""
    remove_outliers(data; nsigma=5) → PhaseData or FrequencyData

Copy of a record with the outliers found by [`find_outliers`](@ref)
replaced by filled values from [`fillgaps`](@ref). A record without
outliers is returned as it is.

An isolated outlier takes the mean of its two neighbors. A run of them is
filled by reflection, with the reservations stated for [`fillgaps`](@ref).

A `PhaseData` record is cleaned as frequency: it is differenced, the
outlying frequency samples are replaced, and the result is summed back to
phase from the first sample. A phase jump is one outlying frequency sample,
so this takes the jump out, where replacing the phase sample it lands on
would leave it in.

Missing samples are filled with the outliers, so a record that has
outliers comes back without gaps. One that has none comes back as it is,
with its gaps.

[SP1065 §10.11](@cite riley-2008-sp1065) and
[Stable32](@cite hamilton-2008-stable32-manual) replace an outlier by a gap,
which the estimators then skip. For a frequency record that is
`y[find_outliers(fd)] .= NaN` on a copy of its samples.

# Examples

```julia
fd = read_frequency("counter_log.csv")
clean = remove_outliers(fd)                # limit of 5
strict = remove_outliers(fd; nsigma=3)
```
"""
function remove_outliers(data::PhaseData; nsigma::Real=5)
    idx = find_outliers(data; nsigma=nsigma)
    isempty(idx) && return data
    # Work in frequency space: a phase step is one outlying frequency sample,
    # so fix it there and integrate back, anchored at the first phase sample.
    pd = any(isnan, data.x) ? fillgaps(data) : data
    x  = pd.x
    T  = eltype(x)
    t0 = pd.tau0
    y  = T[(x[i + 1] - x[i]) / t0 for i in 1:length(x) - 1]
    fd = remove_outliers(FrequencyData(y, t0; source=pd.source); nsigma=nsigma)
    out = Vector{T}(undef, length(x))
    out[1] = x[1]
    @inbounds for i in 1:length(x) - 1
        out[i + 1] = out[i] + fd.y[i] * t0
    end
    return PhaseData(out, t0; source=pd.source, start=pd.start)
end

function remove_outliers(data::FrequencyData; nsigma::Real=5)
    idx = find_outliers(data; nsigma=nsigma)
    isempty(idx) && return data
    y = Vector{Float64}(data.y)
    y[idx] .= NaN
    return fillgaps(FrequencyData(y, data.tau0; source=data.source, start=data.start))
end
