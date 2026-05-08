# Total-family `detrend` kwarg — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` (recommended) or `superpowers:executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `detrend::Symbol` kwarg to all four total-family kernels (`totdev`, `mtotdev`, `htotdev`, `mhtotdev`) and their API wrappers, supporting four named recipes (`:howe`, `:greenhall`, `:linear`, `:legacy`) with canonical per-kernel defaults.

**Architecture:** Per-recipe internal functions (`_<kernel>_<recipe>`) invoked by a per-kernel dispatcher (`_<kernel>_core`). Dispatchers validate the recipe and route to the named helper. Aliases (e.g. `:legacy` MTOT → `_mtotdev_greenhall`) keep code paths minimal. API wrappers pass `detrend` through to the core dispatcher.

**Tech Stack:** Julia 1.11 (`SigmaTauStability` subpackage). No new dependencies. Tests use `Test`, `Random`, `FFTW` (already in test extras).

**Spec:** [`docs/superpowers/specs/2026-05-07-detrend-kwarg-design.md`](../specs/2026-05-07-detrend-kwarg-design.md)

**PR-split suggestion** (set explicitly at execution time):

- PR-A: Phase 1 (infrastructure, no behavior change). All 339 existing tests pass unchanged.
- PR-B: Phase 2 + 3 (new recipes, smoke-tested, defaults still `:legacy`-compatible). Adds new tests, doesn't break old ones.
- PR-C: Phase 4 (default switch — breaking change). Tightens Stable32 TOTDEV cross-val. Updates legacy_kernels.jl test calls. CHANGELOG.

---

## File structure

| File | Action | Responsibility |
|---|---|---|
| `lib/SigmaTauStability/src/core/total.jl` | Modify (refactor + new recipes) | Per-kernel dispatchers + per-recipe helpers |
| `lib/SigmaTauStability/src/api/total.jl` | Modify | Pass `detrend` kwarg through to core |
| `lib/SigmaTauStability/test/legacy_kernels.jl` | No change | Verbatim MATLAB-era reference; unchanged |
| `lib/SigmaTauStability/test/runtests.jl` | Modify | Update legacy parity calls; tighten Stable32 TOT; add cross-recipe + smoke testsets |
| `CHANGELOG.md` | Modify | `[Unreleased] → Changed` entry for the breaking change |
| `TODO.md` | Modify | New entries for MHTOT EDF Monte Carlo + allantools TOT tightening follow-up |

---

## Phase 1 — Infrastructure (refactor + kwarg, no behavior change)

Goal: each `_*_core` becomes a thin dispatcher that delegates to a named `_*_<recipe>` helper. Default recipe = current behavior. All 339 existing tests pass with zero modification.

### Task 1.1: TOTDEV refactor

**Files:**
- Modify: `lib/SigmaTauStability/src/core/total.jl:8-65` (current `_totdev_core` body)

- [ ] **Step 1: Read the existing `_totdev_core` end-to-end**

Open `lib/SigmaTauStability/src/core/total.jl` and confirm the current body (lines 8–65) is the global-LS-detrend + endpoint-mean-flip-reflect implementation. Familiarize yourself with the existing algorithm before refactoring.

- [ ] **Step 2: Rename the existing kernel body to `_totdev_legacy`**

Replace the existing `_totdev_core` definition (lines 3–65) with:

```julia
"""
    _totdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:legacy) → Vector{Float64}

Computes the Total Deviation (TOTDEV) for a set of averaging factors `m`.

`detrend` selects the boundary-handling recipe:
- `:howe` — no detrend, mean-flip endpoint reflection (Howe 1995, NIST SP1065 eqn 25)
- `:greenhall` — per-window half-mean slope removal + time-reverse extension (Greenhall 2003)
- `:linear` — per-window full LS detrend + time-reverse extension
- `:legacy` — current SigmaTau behavior: global LS detrend + mean-flip reflection
"""
function _totdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:legacy)
    detrend === :legacy && return _totdev_legacy(x, m_values, tau0)
    throw(ArgumentError("unknown detrend recipe: $detrend; valid: :howe, :greenhall, :linear, :legacy"))
end

function _totdev_legacy(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    # Linear detrend of the whole vector (analytic LS sums)
    N_float = Float64(N)
    sum_i = (N_float * (N_float + 1.0)) / 2.0
    sum_i2 = (N_float * (N_float + 1.0) * (2.0*N_float + 1.0)) / 6.0
    delta = N_float * sum_i2 - sum_i^2

    sum_x = sum(x)
    sum_ix = 0.0
    @inbounds @simd for i in 1:N
        sum_ix += i * x[i]
    end

    a = (sum_x * sum_i2 - sum_ix * sum_i) / delta
    b = (N_float * sum_ix - sum_x * sum_i) / delta

    xd = Vector{Float64}(undef, N)
    @inbounds @simd for i in 1:N
        xd[i] = x[i] - (a + b * i)
    end

    # Mean-flip endpoint reflection: x_star of length 3N-4
    x_star = Vector{Float64}(undef, 3N - 4)
    @inbounds for i in 1:N-2
        x_star[i] = 2.0*xd[1] - xd[i+1]
    end
    @inbounds for i in 1:N
        x_star[N-2+i] = xd[i]
    end
    @inbounds for i in 1:N-2
        x_star[2N-2+i] = 2.0*xd[N] - xd[N-i]
    end

    off = N - 2
    for (k, m) in enumerate(m_values)
        D = 0.0
        count = 0
        @inbounds @simd for i in 1:N
            lo = off + i
            hi = off + i + 2m
            if hi <= length(x_star)
                d2 = x_star[hi] - 2.0*x_star[off + i + m] + x_star[lo]
                D += d2^2
                count += 1
            end
        end

        if count == 0
            devs[k] = NaN
        else
            devs[k] = sqrt(D / (2.0 * (N - 2) * Float64(m)^2 * tau0^2))
        end
    end

    return devs
end
```

- [ ] **Step 3: Run the existing test suite**

```
julia> ]
pkg> test SigmaTauStability
```

Expected: `339/339` pass (or whatever your current baseline is). The dispatcher routes `:legacy` (default) to the renamed helper; behavior is identical.

- [ ] **Step 4: Commit**

```
git add lib/SigmaTauStability/src/core/total.jl
git commit -m "refactor: extract _totdev_legacy from _totdev_core, add detrend dispatcher"
```

### Task 1.2: MTOTDEV refactor

**Files:**
- Modify: `lib/SigmaTauStability/src/core/total.jl:67-140` (current `_mtotdev_core` body)

- [ ] **Step 1: Rename the existing kernel body to `_mtotdev_greenhall`**

The current `_mtotdev_core` body IS the canonical Greenhall recipe (per-window half-mean slope removal + time-reverse extension + modified 2nd-diff). Rename it.

Replace the existing `_mtotdev_core` definition with:

```julia
"""
    _mtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:greenhall) → Vector{Float64}

Computes the Modified Total Deviation (MTOTDEV).

`detrend` selects the boundary-handling recipe (see `_totdev_core` docstring).
For MTOTDEV, `:legacy` is an alias for `:greenhall` (current implementation).
"""
function _mtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:greenhall)
    if detrend === :greenhall || detrend === :legacy
        return _mtotdev_greenhall(x, m_values, tau0)
    end
    throw(ArgumentError("unknown detrend recipe: $detrend; valid: :howe, :greenhall, :linear, :legacy"))
end

function _mtotdev_greenhall(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    # [Move the existing _mtotdev_core body here verbatim — it already implements
    #  per-window half-mean slope removal + time-reverse extension + modified 2nd-diff,
    #  which is the Greenhall recipe.]
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    max_m = isempty(m_values) ? 0 : maximum(m_values)
    max_seg = 3 * max_m
    ext = Vector{Float64}(undef, 3 * max_seg)
    cs = Vector{Float64}(undef, 3 * max_seg + 1)

    for (k, m) in enumerate(m_values)
        nsubs = N - 3m + 1
        if nsubs < 1
            devs[k] = NaN
            continue
        end

        seg_len = 3m
        outer_sum = 0.0

        for n in 1:nsubs
            half_n = seg_len / 2.0
            if m == 1
                slope = (x[n+2] - x[n]) / (2.0 * tau0)
            else
                hi = floor(Int, half_n)
                s1 = 0.0
                @inbounds @simd for i in 1:hi
                    s1 += x[n-1+i]
                end
                s1 /= hi

                s2 = 0.0
                @inbounds @simd for i in hi+1:seg_len
                    s2 += x[n-1+i]
                end
                s2 /= (seg_len - hi)
                slope = (s2 - s1) / (half_n * tau0)
            end

            @inbounds for j in 1:seg_len
                val = x[n-1+j] - slope * tau0 * (j - 1)
                rev_val = x[n-1 + seg_len - j + 1] - slope * tau0 * (seg_len - j)

                ext[j] = rev_val
                ext[seg_len + j] = val
                ext[2seg_len + j] = rev_val
            end

            cs[1] = 0.0
            @inbounds for j in 1:3seg_len
                cs[j+1] = cs[j] + ext[j]
            end

            block_sum = 0.0
            @inbounds @simd for j in 0:(6m - 1)
                a1 = (cs[j+m+1]  - cs[j+1])
                a2 = (cs[j+2m+1] - cs[j+m+1])
                a3 = (cs[j+3m+1] - cs[j+2m+1])
                d2 = (a3 - 2.0*a2 + a1) / m
                block_sum += d2^2
            end
            outer_sum += block_sum / (6.0 * m)
        end

        devs[k] = sqrt(outer_sum / (2.0 * Float64(m)^2 * tau0^2 * nsubs))
    end

    return devs
end
```

- [ ] **Step 2: Run tests, expect pass**

```
pkg> test SigmaTauStability
```

339/339 pass. `:greenhall` (default) and `:legacy` both route to the same helper.

- [ ] **Step 3: Commit**

```
git add lib/SigmaTauStability/src/core/total.jl
git commit -m "refactor: extract _mtotdev_greenhall, alias :legacy to it"
```

### Task 1.3: HTOTDEV refactor

**Files:**
- Modify: `lib/SigmaTauStability/src/core/total.jl:142-234`

- [ ] **Step 1: Rename HTOTDEV body to `_htotdev_greenhall`**

The current `_htotdev_core` body is the canonical Greenhall recipe on the frequency series (`y = diff(x)/tau0`), with per-window half-mean slope + time-reverse extension + third-difference. Rename and add dispatcher.

Replace the existing `_htotdev_core` definition with:

```julia
"""
    _htotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:greenhall) → Vector{Float64}

Computes the Hadamard Total Deviation (HTOTDEV).

`detrend` selects the boundary-handling recipe (see `_totdev_core`). For HTOTDEV,
`:legacy` is an alias for `:greenhall`. The recipe operates on the frequency
series `y = diff(x) / tau0`.
"""
function _htotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:greenhall)
    if detrend === :greenhall || detrend === :legacy
        return _htotdev_greenhall(x, m_values, tau0)
    end
    throw(ArgumentError("unknown detrend recipe: $detrend; valid: :howe, :greenhall, :linear, :legacy"))
end

function _htotdev_greenhall(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    # [Move the existing _htotdev_core body here verbatim.]
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    y = Vector{Float64}(undef, N-1)
    @inbounds @simd for i in 1:N-1
        y[i] = (x[i+1] - x[i]) / tau0
    end
    Ny = length(y)

    max_m = isempty(m_values) ? 0 : maximum(m_values)
    max_seg = 3 * max_m
    ext = Vector{Float64}(undef, 3 * max_seg)
    cs = Vector{Float64}(undef, 3 * max_seg + 1)

    for (k, m) in enumerate(m_values)
        if m == 1
            L = N - 3
            if L <= 0
                devs[k] = NaN
                continue
            end
            sum_sq = 0.0
            @inbounds @simd for i in 1:L
                d3 = x[i+3] - 3.0*x[i+2] + 3.0*x[i+1] - x[i]
                sum_sq += d3^2
            end
            devs[k] = sqrt(sum_sq / (6.0 * L * tau0^2))
            continue
        end

        n_iter = Ny - 3m + 1
        if n_iter < 1
            devs[k] = NaN
            continue
        end

        seg_len = 3m
        dev_sum = 0.0

        for i in 0:(n_iter - 1)
            hi = floor(Int, seg_len / 2)
            lo_start = ceil(Int, seg_len / 2) + 1

            s1 = 0.0
            @inbounds @simd for j in 1:hi
                s1 += y[i+j]
            end
            m1 = s1 / hi

            s2 = 0.0
            @inbounds @simd for j in lo_start:seg_len
                s2 += y[i+j]
            end
            m2 = s2 / (seg_len - lo_start + 1)

            slope = isodd(seg_len) ? (m2 - m1) / (0.5*(seg_len - 1) + 1.0) : (m2 - m1) / (0.5*seg_len)
            mid = floor(seg_len / 2)

            @inbounds for j in 1:seg_len
                val = y[i+j] - slope * (j - 1 - mid)
                rev_val = y[i + seg_len - j + 1] - slope * (seg_len - j - mid)

                ext[j] = rev_val
                ext[seg_len + j] = val
                ext[2seg_len + j] = rev_val
            end

            cs[1] = 0.0
            @inbounds for j in 1:3seg_len
                cs[j+1] = cs[j] + ext[j]
            end

            sq = 0.0
            @inbounds @simd for j in 0:(6m - 1)
                h1 = (cs[j+m+1]  - cs[j+1])
                h2 = (cs[j+2m+1] - cs[j+m+1])
                h3 = (cs[j+3m+1] - cs[j+2m+1])
                sq += ((h3 - 2.0*h2 + h1) / m)^2
            end
            dev_sum += sq / (6.0 * m)
        end

        devs[k] = sqrt(dev_sum / (6.0 * n_iter))
    end

    return devs
end
```

- [ ] **Step 2: Run tests, expect 339/339 pass**

```
pkg> test SigmaTauStability
```

- [ ] **Step 3: Commit**

```
git add lib/SigmaTauStability/src/core/total.jl
git commit -m "refactor: extract _htotdev_greenhall, alias :legacy to it"
```

### Task 1.4: MHTOTDEV refactor

**Files:**
- Modify: `lib/SigmaTauStability/src/core/total.jl:236-321`

- [ ] **Step 1: Rename MHTOTDEV body to `_mhtotdev_linear`**

Current MHTOTDEV body = per-window full LS detrend + time-reverse extension + averaged third-difference. That IS the `:linear` recipe shape. Rename and add dispatcher.

Default for MHTOTDEV in this phase is `:linear` (preserves current behavior). The default switch to `:greenhall` happens in Phase 4.

Replace the existing `_mhtotdev_core` definition with:

```julia
"""
    _mhtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:linear) → Vector{Float64}

Computes the Modified Hadamard Total Deviation (MHTOTDEV).

`detrend` selects the boundary-handling recipe (see `_totdev_core`). MHTOTDEV
is novel to SigmaTau; `:legacy` aliases to `:linear` (current implementation).
The default is `:linear` in this phase; switches to `:greenhall` in Phase 4.
"""
function _mhtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:linear)
    if detrend === :linear || detrend === :legacy
        return _mhtotdev_linear(x, m_values, tau0)
    end
    throw(ArgumentError("unknown detrend recipe: $detrend; valid: :howe, :greenhall, :linear, :legacy"))
end

function _mhtotdev_linear(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    # [Move the existing _mhtotdev_core body here verbatim.]
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    max_m = isempty(m_values) ? 0 : maximum(m_values)
    max_Lp = 3 * max_m + 1
    ext_len = 3 * max_Lp
    L3_max = ext_len - 3 * max_m
    ext = Vector{Float64}(undef, ext_len)
    d3_vec = Vector{Float64}(undef, L3_max)
    S = Vector{Float64}(undef, L3_max + 1)

    for (k, m) in enumerate(m_values)
        if m < 1
            devs[k] = NaN
            continue
        end
        nsubs = N - 4m + 1
        if nsubs < 1
            devs[k] = NaN
            continue
        end

        Lp = 3m + 1
        L3 = 3Lp - 3m

        total_sum = 0.0
        for n in 1:nsubs
            Lp_float = Float64(Lp)
            sum_i = (Lp_float * (Lp_float + 1.0)) / 2.0
            sum_i2 = (Lp_float * (Lp_float + 1.0) * (2.0*Lp_float + 1.0)) / 6.0
            delta = Lp_float * sum_i2 - sum_i^2

            sum_x = 0.0
            sum_ix = 0.0
            @inbounds @simd for j in 1:Lp
                val = x[n-1+j]
                sum_x += val
                sum_ix += j * val
            end

            a = (sum_x * sum_i2 - sum_ix * sum_i) / delta
            b = (Lp_float * sum_ix - sum_x * sum_i) / delta

            @inbounds for j in 1:Lp
                val = x[n-1+j] - (a + b * j)
                rev_val = x[n-1 + Lp - j + 1] - (a + b * (Lp - j + 1))

                ext[j] = rev_val
                ext[Lp + j] = val
                ext[2Lp + j] = rev_val
            end

            @inbounds for j in 1:L3
                d3_vec[j] = ext[j] - 3.0*ext[j+m] + 3.0*ext[j+2m] - ext[j+3m]
            end

            S[1] = 0.0
            @inbounds for j in 1:L3
                S[j+1] = S[j] + d3_vec[j]
            end

            n_avg = L3 + 1 - m
            if n_avg > 0
                block_var = 0.0
                @inbounds @simd for j in 1:n_avg
                    block_var += (S[j+m] - S[j])^2
                end
                block_var /= (n_avg * 6.0 * Float64(m)^2)
            else
                block_var = 0.0
            end

            total_sum += block_var
        end

        devs[k] = sqrt(total_sum / (nsubs * Float64(m)^2 * tau0^2))
    end

    return devs
end
```

- [ ] **Step 2: Run tests, expect 339/339 pass**

```
pkg> test SigmaTauStability
```

- [ ] **Step 3: Commit**

```
git add lib/SigmaTauStability/src/core/total.jl
git commit -m "refactor: extract _mhtotdev_linear, alias :legacy to it"
```

### Task 1.5: API wrapper plumbing

**Files:**
- Modify: `lib/SigmaTauStability/src/api/total.jl:8-96` (all four wrappers)

- [ ] **Step 1: Add `detrend` kwarg to each wrapper**

For each of `totdev`, `mtotdev`, `htotdev`, `mhtotdev`, add a `detrend::Symbol` kwarg with a default that matches the corresponding `_*_core` default. Pass through to the core call.

Example for `totdev`:

```julia
"""
    totdev(data::PhaseData, m_values::Vector{Int}; detrend::Symbol=:legacy, calc_ci::Bool=true, confidence::Float64=0.95)

Computes the Total Deviation for the given PhaseData. See `_totdev_core` for
the meaning of `detrend`.
"""
function totdev(data::PhaseData, m_values::Vector{Int}; detrend::Symbol=:legacy, calc_ci::Bool=true, confidence::Float64=0.95)
    raw_devs = _totdev_core(data.x, m_values, data.tau0; detrend=detrend)
    taus = m_values .* data.tau0
    T = (length(data.x) - 1) * data.tau0

    if !calc_ci
        return StabilityResult(:totdev, taus, raw_devs, Symbol[], Float64[], Float64[], Float64[])
    end

    noises = identify_noise(data.x, m_values, dmin=0, dmax=2)
    B = bias_correction(noises, :totvar, taus, T)
    biased_devs = raw_devs ./ B

    edfs = calculate_edf(:totdev, biased_devs, noises, m_values, taus, length(data.x), T)
    lower, upper = confidence_intervals(biased_devs, edfs, noises, length(data.x), confidence)

    return StabilityResult(:totdev, taus, biased_devs, noises, lower, upper, edfs)
end
```

Apply the same pattern to `mtotdev` (default `:greenhall`), `htotdev` (default `:greenhall`), and `mhtotdev` (default `:linear`).

Also update the `FrequencyData` entry points at the bottom of `api/total.jl` to pass `detrend` through:

```julia
totdev(data::FrequencyData, m_values::Vector{Int}; kwargs...)   = totdev(_freq_to_phase(data),   m_values; kwargs...)
mtotdev(data::FrequencyData, m_values::Vector{Int}; kwargs...)  = mtotdev(_freq_to_phase(data),  m_values; kwargs...)
htotdev(data::FrequencyData, m_values::Vector{Int}; kwargs...)  = htotdev(_freq_to_phase(data),  m_values; kwargs...)
mhtotdev(data::FrequencyData, m_values::Vector{Int}; kwargs...) = mhtotdev(_freq_to_phase(data), m_values; kwargs...)
```

`kwargs...` captures and forwards the `detrend` kwarg automatically — no change needed to those lines.

- [ ] **Step 2: Run tests, expect 339/339 pass**

```
pkg> test SigmaTauStability
```

- [ ] **Step 3: Commit**

```
git add lib/SigmaTauStability/src/api/total.jl
git commit -m "feat: detrend kwarg on totdev/mtotdev/htotdev/mhtotdev API wrappers"
```

**Phase 1 complete.** Push the branch:

```
git push origin track-a4-detrend
```

This is a natural PR-A boundary if splitting.

---

## Phase 2 — TOTDEV recipes

Add `:howe`, `:greenhall`, `:linear` recipes to TOTDEV. Each is TDD: failing test → implement → pass → commit.

### Task 2.1: TOTDEV `:howe` (canonical SP1065 eqn 25)

**Files:**
- Modify: `lib/SigmaTauStability/src/core/total.jl` (add `_totdev_howe` + dispatcher branch)
- Modify: `lib/SigmaTauStability/test/runtests.jl` (new testset)

- [ ] **Step 1: Write the failing test**

In `lib/SigmaTauStability/test/runtests.jl`, after the existing `@testset "Stable32 cross-validation"` block (around line 300), add a new testset:

```julia
@testset "TOTDEV :howe matches Stable32 tightly" begin
    # SP1065 eqn 25 reference: no detrend, mean-flip endpoint reflection.
    # Should match Stable32's TOTDEV output at rtol=1e-4 (vs the rtol=0.15
    # boundary-policy floor seen with the :legacy global-LS detrend).
    ref_dir = joinpath(@__DIR__, "..", "..", "..", "reference", "validation")
    dat_path = joinpath(ref_dir, "stable32gen.DAT")
    csv_path = joinpath(ref_dir, "stable32out", "stable32_data_full.csv")

    if !isfile(dat_path) || !isfile(csv_path)
        @warn "Stable32 fixtures not present, skipping :howe TOTDEV tightness test"
    else
        lines = readlines(dat_path)
        x = parse.(Float64, strip.(lines[11:end]))
        @test length(x) == 8192
        tau0 = 1.0

        rows = [split(line, ',') for line in readlines(csv_path)[2:end]]
        n_checked = 0
        for row in rows
            length(row) < 7 && continue
            row[1] == "Total" || continue
            m = parse(Int, row[2])
            sigma_ref = parse(Float64, row[7])

            got = SigmaTauStability._totdev_core(x, [m], tau0; detrend=:howe)[1]
            @test got ≈ sigma_ref rtol=1e-4
            n_checked += 1
        end
        @test n_checked >= 5
    end
end
```

- [ ] **Step 2: Run the test, expect fail**

```
pkg> test SigmaTauStability
```

Expected: ArgumentError "unknown detrend recipe: howe" — dispatcher rejects `:howe` until we implement it.

- [ ] **Step 3: Implement `_totdev_howe`**

In `lib/SigmaTauStability/src/core/total.jl`, after `_totdev_legacy`, add:

```julia
function _totdev_howe(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    # SP1065 eqn 25: no detrend, mean-flip endpoint reflection.
    # x_star of length 3N-4: [reflected_left | x | reflected_right].
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    x_star = Vector{Float64}(undef, 3N - 4)
    @inbounds for i in 1:N-2
        x_star[i] = 2.0*x[1] - x[i+1]
    end
    @inbounds for i in 1:N
        x_star[N-2+i] = x[i]
    end
    @inbounds for i in 1:N-2
        x_star[2N-2+i] = 2.0*x[N] - x[N-i]
    end

    off = N - 2
    for (k, m) in enumerate(m_values)
        D = 0.0
        count = 0
        @inbounds @simd for i in 1:N
            lo = off + i
            hi = off + i + 2m
            if hi <= length(x_star)
                d2 = x_star[hi] - 2.0*x_star[off + i + m] + x_star[lo]
                D += d2^2
                count += 1
            end
        end

        devs[k] = count == 0 ? NaN : sqrt(D / (2.0 * (N - 2) * Float64(m)^2 * tau0^2))
    end

    return devs
end
```

Then update the `_totdev_core` dispatcher to add the `:howe` branch:

```julia
function _totdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:legacy)
    detrend === :legacy && return _totdev_legacy(x, m_values, tau0)
    detrend === :howe && return _totdev_howe(x, m_values, tau0)
    throw(ArgumentError("unknown detrend recipe: $detrend; valid: :howe, :greenhall, :linear, :legacy"))
end
```

- [ ] **Step 4: Run the test, expect pass**

```
pkg> test SigmaTauStability
```

The new testset passes; existing 339 + new ones all green.

- [ ] **Step 5: Commit**

```
git add lib/SigmaTauStability/src/core/total.jl lib/SigmaTauStability/test/runtests.jl
git commit -m "feat: TOTDEV :howe recipe (SP1065 eqn 25, matches Stable32 at rtol=1e-4)"
```

### Task 2.2: TOTDEV `:greenhall` (per-window half-mean + time-reverse)

**Files:**
- Modify: `lib/SigmaTauStability/src/core/total.jl`
- Modify: `lib/SigmaTauStability/test/runtests.jl`

- [ ] **Step 1: Write the failing smoke test**

In `runtests.jl`, after the previous testset, add:

```julia
@testset "TOTDEV :greenhall smoke" begin
    using Random
    Random.seed!(20260508)
    N    = 1024
    tau0 = 1.0
    ms   = [1, 2, 4, 8, 16]
    # WHFM (alpha = 0) — mid-spectrum noise type
    x = _gen_powerlaw_phase(0.0, N; tau0=tau0)

    devs = SigmaTauStability._totdev_core(x, ms, tau0; detrend=:greenhall)
    @test length(devs) == length(ms)
    @test all(isfinite, devs)
    @test all(>(0), devs)

    # Should be within an order of magnitude of :legacy on the same fixture.
    devs_legacy = SigmaTauStability._totdev_core(x, ms, tau0; detrend=:legacy)
    @test all(0.1 .<= devs ./ devs_legacy .<= 10.0)
end
```

- [ ] **Step 2: Run the test, expect fail (unknown recipe)**

```
pkg> test SigmaTauStability
```

- [ ] **Step 3: Implement `_totdev_greenhall`**

In `core/total.jl`, after `_totdev_howe`, add:

```julia
function _totdev_greenhall(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    # Per-window half-mean slope removal + time-reverse extension + standard 2nd-diff.
    # Per Greenhall 2003. Uses a 3m-sample window per (n, m) pair.
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    max_m = isempty(m_values) ? 0 : maximum(m_values)
    max_seg = 3 * max_m
    ext = Vector{Float64}(undef, 3 * max_seg)

    for (k, m) in enumerate(m_values)
        nsubs = N - 3m + 1
        if nsubs < 1
            devs[k] = NaN
            continue
        end

        seg_len = 3m
        outer_sum = 0.0

        for n in 1:nsubs
            half_n = seg_len / 2.0
            if m == 1
                slope = (x[n+2] - x[n]) / (2.0 * tau0)
            else
                hi = floor(Int, half_n)
                s1 = 0.0
                @inbounds @simd for i in 1:hi
                    s1 += x[n-1+i]
                end
                s1 /= hi

                s2 = 0.0
                @inbounds @simd for i in hi+1:seg_len
                    s2 += x[n-1+i]
                end
                s2 /= (seg_len - hi)
                slope = (s2 - s1) / (half_n * tau0)
            end

            # Detrend window + time-reverse extension
            @inbounds for j in 1:seg_len
                val = x[n-1+j] - slope * tau0 * (j - 1)
                rev_val = x[n-1 + seg_len - j + 1] - slope * tau0 * (seg_len - j)

                ext[j] = rev_val
                ext[seg_len + j] = val
                ext[2seg_len + j] = rev_val
            end

            # Standard 2nd-difference operator (TOTDEV form, not modified) on the
            # extended series. Slide over the original-data positions only.
            block_sum = 0.0
            count = 0
            @inbounds @simd for i in 1:seg_len
                lo = seg_len + i
                hi = seg_len + i + 2m
                if hi <= 3seg_len
                    d2 = ext[hi] - 2.0*ext[lo + m] + ext[lo]
                    block_sum += d2^2
                    count += 1
                end
            end

            outer_sum += count > 0 ? block_sum / (2.0 * count * Float64(m)^2 * tau0^2) : 0.0
        end

        devs[k] = sqrt(outer_sum / nsubs)
    end

    return devs
end
```

Update dispatcher:

```julia
function _totdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:legacy)
    detrend === :legacy && return _totdev_legacy(x, m_values, tau0)
    detrend === :howe && return _totdev_howe(x, m_values, tau0)
    detrend === :greenhall && return _totdev_greenhall(x, m_values, tau0)
    throw(ArgumentError("unknown detrend recipe: $detrend; valid: :howe, :greenhall, :linear, :legacy"))
end
```

- [ ] **Step 4: Run the test, expect pass**

```
pkg> test SigmaTauStability
```

- [ ] **Step 5: Commit**

```
git add lib/SigmaTauStability/src/core/total.jl lib/SigmaTauStability/test/runtests.jl
git commit -m "feat: TOTDEV :greenhall recipe (per-window half-mean + time-reverse)"
```

### Task 2.3: TOTDEV `:linear` (per-window full LS + time-reverse)

**Files:**
- Modify: `lib/SigmaTauStability/src/core/total.jl`
- Modify: `lib/SigmaTauStability/test/runtests.jl`

- [ ] **Step 1: Write the failing smoke test**

In `runtests.jl`, after the previous testset:

```julia
@testset "TOTDEV :linear smoke" begin
    using Random
    Random.seed!(20260508)
    N    = 1024
    tau0 = 1.0
    ms   = [1, 2, 4, 8, 16]
    x = _gen_powerlaw_phase(0.0, N; tau0=tau0)

    devs = SigmaTauStability._totdev_core(x, ms, tau0; detrend=:linear)
    @test length(devs) == length(ms)
    @test all(isfinite, devs)
    @test all(>(0), devs)
    devs_legacy = SigmaTauStability._totdev_core(x, ms, tau0; detrend=:legacy)
    @test all(0.1 .<= devs ./ devs_legacy .<= 10.0)
end
```

- [ ] **Step 2: Run, expect fail**

```
pkg> test SigmaTauStability
```

- [ ] **Step 3: Implement `_totdev_linear`**

In `core/total.jl`, after `_totdev_greenhall`, add:

```julia
function _totdev_linear(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    # Per-window full LS detrend + time-reverse extension + standard 2nd-diff.
    # Same structure as _totdev_greenhall but uses analytic-LS slope+intercept
    # instead of half-mean slope.
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    max_m = isempty(m_values) ? 0 : maximum(m_values)
    max_seg = 3 * max_m
    ext = Vector{Float64}(undef, 3 * max_seg)

    for (k, m) in enumerate(m_values)
        nsubs = N - 3m + 1
        if nsubs < 1
            devs[k] = NaN
            continue
        end

        seg_len = 3m
        L_float = Float64(seg_len)
        sum_i = (L_float * (L_float + 1.0)) / 2.0
        sum_i2 = (L_float * (L_float + 1.0) * (2.0*L_float + 1.0)) / 6.0
        delta = L_float * sum_i2 - sum_i^2

        outer_sum = 0.0

        for n in 1:nsubs
            sum_x = 0.0
            sum_ix = 0.0
            @inbounds @simd for j in 1:seg_len
                v = x[n-1+j]
                sum_x += v
                sum_ix += j * v
            end

            a = (sum_x * sum_i2 - sum_ix * sum_i) / delta
            b = (L_float * sum_ix - sum_x * sum_i) / delta

            @inbounds for j in 1:seg_len
                val = x[n-1+j] - (a + b * j)
                rev_val = x[n-1 + seg_len - j + 1] - (a + b * (seg_len - j + 1))

                ext[j] = rev_val
                ext[seg_len + j] = val
                ext[2seg_len + j] = rev_val
            end

            block_sum = 0.0
            count = 0
            @inbounds @simd for i in 1:seg_len
                lo = seg_len + i
                hi = seg_len + i + 2m
                if hi <= 3seg_len
                    d2 = ext[hi] - 2.0*ext[lo + m] + ext[lo]
                    block_sum += d2^2
                    count += 1
                end
            end

            outer_sum += count > 0 ? block_sum / (2.0 * count * Float64(m)^2 * tau0^2) : 0.0
        end

        devs[k] = sqrt(outer_sum / nsubs)
    end

    return devs
end
```

Update dispatcher:

```julia
function _totdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:legacy)
    detrend === :legacy && return _totdev_legacy(x, m_values, tau0)
    detrend === :howe && return _totdev_howe(x, m_values, tau0)
    detrend === :greenhall && return _totdev_greenhall(x, m_values, tau0)
    detrend === :linear && return _totdev_linear(x, m_values, tau0)
    throw(ArgumentError("unknown detrend recipe: $detrend; valid: :howe, :greenhall, :linear, :legacy"))
end
```

- [ ] **Step 4: Run, expect pass**

```
pkg> test SigmaTauStability
```

- [ ] **Step 5: Commit**

```
git add lib/SigmaTauStability/src/core/total.jl lib/SigmaTauStability/test/runtests.jl
git commit -m "feat: TOTDEV :linear recipe (per-window LS + time-reverse)"
```

---

## Phase 3 — MTOT / HTOT / MHTOT new recipes

Each kernel needs `:howe` and `:linear` (MTOT/HTOT) or `:howe` and `:greenhall` (MHTOT) added. Same TDD pattern.

### Task 3.1: MTOTDEV `:howe`

**Files:**
- Modify: `lib/SigmaTauStability/src/core/total.jl`
- Modify: `lib/SigmaTauStability/test/runtests.jl`

- [ ] **Step 1: Write smoke test**

```julia
@testset "MTOTDEV :howe smoke" begin
    using Random
    Random.seed!(20260508)
    N    = 1024
    tau0 = 1.0
    ms   = [1, 2, 4, 8, 16]
    x = _gen_powerlaw_phase(0.0, N; tau0=tau0)

    devs = SigmaTauStability._mtotdev_core(x, ms, tau0; detrend=:howe)
    @test length(devs) == length(ms)
    @test all(isfinite, devs)
    @test all(>(0), devs)
    devs_legacy = SigmaTauStability._mtotdev_core(x, ms, tau0; detrend=:legacy)
    @test all(0.1 .<= devs ./ devs_legacy .<= 10.0)
end
```

- [ ] **Step 2: Run, expect fail**

- [ ] **Step 3: Implement `_mtotdev_howe`**

```julia
function _mtotdev_howe(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    # No detrend, global mean-flip endpoint reflection (length 3N-4).
    # Then apply the modified 2nd-diff operator on the extended series.
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    x_star = Vector{Float64}(undef, 3N - 4)
    @inbounds for i in 1:N-2
        x_star[i] = 2.0*x[1] - x[i+1]
    end
    @inbounds for i in 1:N
        x_star[N-2+i] = x[i]
    end
    @inbounds for i in 1:N-2
        x_star[2N-2+i] = 2.0*x[N] - x[N-i]
    end

    L = length(x_star)
    cs = pushfirst!(cumsum(x_star), 0.0)

    for (k, m) in enumerate(m_values)
        nsubs = L - 3m + 1
        if nsubs < 1
            devs[k] = NaN
            continue
        end
        outer_sum = 0.0
        @inbounds @simd for j in 0:(nsubs - 1)
            a1 = cs[j+m+1]   - cs[j+1]
            a2 = cs[j+2m+1]  - cs[j+m+1]
            a3 = cs[j+3m+1]  - cs[j+2m+1]
            d2 = (a3 - 2.0*a2 + a1) / m
            outer_sum += d2^2
        end
        devs[k] = sqrt(outer_sum / (2.0 * Float64(m)^2 * tau0^2 * nsubs))
    end

    return devs
end
```

Update dispatcher:

```julia
function _mtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:greenhall)
    if detrend === :greenhall || detrend === :legacy
        return _mtotdev_greenhall(x, m_values, tau0)
    end
    detrend === :howe && return _mtotdev_howe(x, m_values, tau0)
    throw(ArgumentError("unknown detrend recipe: $detrend; valid: :howe, :greenhall, :linear, :legacy"))
end
```

- [ ] **Step 4: Run, expect pass**

- [ ] **Step 5: Commit**

```
git add lib/SigmaTauStability/src/core/total.jl lib/SigmaTauStability/test/runtests.jl
git commit -m "feat: MTOTDEV :howe recipe (no detrend, mean-flip reflect)"
```

### Task 3.2: MTOTDEV `:linear`

**Files:**
- Modify: `lib/SigmaTauStability/src/core/total.jl`
- Modify: `lib/SigmaTauStability/test/runtests.jl`

- [ ] **Step 1: Write smoke test (mirrors 3.1, swap :howe → :linear)**

```julia
@testset "MTOTDEV :linear smoke" begin
    using Random
    Random.seed!(20260508)
    N    = 1024
    tau0 = 1.0
    ms   = [1, 2, 4, 8, 16]
    x = _gen_powerlaw_phase(0.0, N; tau0=tau0)

    devs = SigmaTauStability._mtotdev_core(x, ms, tau0; detrend=:linear)
    @test length(devs) == length(ms)
    @test all(isfinite, devs)
    @test all(>(0), devs)
    devs_legacy = SigmaTauStability._mtotdev_core(x, ms, tau0; detrend=:legacy)
    @test all(0.1 .<= devs ./ devs_legacy .<= 10.0)
end
```

- [ ] **Step 2: Run, expect fail**

- [ ] **Step 3: Implement `_mtotdev_linear`**

Derive from `_mtotdev_greenhall` by replacing the half-mean slope estimation (lines computing `s1`, `s2`, `slope`) with a closed-form analytic LS fit `(a, b)`, exactly as in `_mhtotdev_linear`. Then the per-window detrend uses `x[j] - (a + b*j)` instead of `x[j] - slope*tau0*(j-1)`. The time-reverse extension and modified 2nd-diff operator are unchanged.

```julia
function _mtotdev_linear(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    # Per-window full LS detrend + time-reverse extension + modified 2nd-diff.
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    max_m = isempty(m_values) ? 0 : maximum(m_values)
    max_seg = 3 * max_m
    ext = Vector{Float64}(undef, 3 * max_seg)
    cs = Vector{Float64}(undef, 3 * max_seg + 1)

    for (k, m) in enumerate(m_values)
        nsubs = N - 3m + 1
        if nsubs < 1
            devs[k] = NaN
            continue
        end

        seg_len = 3m
        L_float = Float64(seg_len)
        sum_i = (L_float * (L_float + 1.0)) / 2.0
        sum_i2 = (L_float * (L_float + 1.0) * (2.0*L_float + 1.0)) / 6.0
        delta = L_float * sum_i2 - sum_i^2

        outer_sum = 0.0

        for n in 1:nsubs
            sum_x = 0.0
            sum_ix = 0.0
            @inbounds @simd for j in 1:seg_len
                v = x[n-1+j]
                sum_x += v
                sum_ix += j * v
            end

            a = (sum_x * sum_i2 - sum_ix * sum_i) / delta
            b = (L_float * sum_ix - sum_x * sum_i) / delta

            @inbounds for j in 1:seg_len
                val = x[n-1+j] - (a + b * j)
                rev_val = x[n-1 + seg_len - j + 1] - (a + b * (seg_len - j + 1))

                ext[j] = rev_val
                ext[seg_len + j] = val
                ext[2seg_len + j] = rev_val
            end

            cs[1] = 0.0
            @inbounds for j in 1:3seg_len
                cs[j+1] = cs[j] + ext[j]
            end

            block_sum = 0.0
            @inbounds @simd for j in 0:(6m - 1)
                a1 = (cs[j+m+1]  - cs[j+1])
                a2 = (cs[j+2m+1] - cs[j+m+1])
                a3 = (cs[j+3m+1] - cs[j+2m+1])
                d2 = (a3 - 2.0*a2 + a1) / m
                block_sum += d2^2
            end
            outer_sum += block_sum / (6.0 * m)
        end

        devs[k] = sqrt(outer_sum / (2.0 * Float64(m)^2 * tau0^2 * nsubs))
    end

    return devs
end
```

Update dispatcher to add `:linear` branch:

```julia
function _mtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:greenhall)
    if detrend === :greenhall || detrend === :legacy
        return _mtotdev_greenhall(x, m_values, tau0)
    end
    detrend === :howe && return _mtotdev_howe(x, m_values, tau0)
    detrend === :linear && return _mtotdev_linear(x, m_values, tau0)
    throw(ArgumentError("unknown detrend recipe: $detrend; valid: :howe, :greenhall, :linear, :legacy"))
end
```

- [ ] **Step 4: Run, expect pass**

- [ ] **Step 5: Commit**

```
git commit -m "feat: MTOTDEV :linear recipe (per-window LS + time-reverse)"
```

### Task 3.3: HTOTDEV `:howe`

**Files:**
- Modify: `lib/SigmaTauStability/src/core/total.jl`
- Modify: `lib/SigmaTauStability/test/runtests.jl`

- [ ] **Step 1: Write smoke test (mirrors 3.1 for HTOT)**

```julia
@testset "HTOTDEV :howe smoke" begin
    using Random
    Random.seed!(20260508)
    N    = 1024
    tau0 = 1.0
    ms   = [1, 2, 4, 8, 16]
    x = _gen_powerlaw_phase(0.0, N; tau0=tau0)

    devs = SigmaTauStability._htotdev_core(x, ms, tau0; detrend=:howe)
    @test length(devs) == length(ms)
    @test all(isfinite, devs)
    @test all(>(0), devs)
    devs_legacy = SigmaTauStability._htotdev_core(x, ms, tau0; detrend=:legacy)
    @test all(0.1 .<= devs ./ devs_legacy .<= 10.0)
end
```

- [ ] **Step 2: Run, expect fail**

- [ ] **Step 3: Implement `_htotdev_howe`**

```julia
function _htotdev_howe(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    # No detrend, global mean-flip reflection on phase, then convert to
    # frequency series and apply third-difference operator.
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    x_star = Vector{Float64}(undef, 3N - 4)
    @inbounds for i in 1:N-2
        x_star[i] = 2.0*x[1] - x[i+1]
    end
    @inbounds for i in 1:N
        x_star[N-2+i] = x[i]
    end
    @inbounds for i in 1:N-2
        x_star[2N-2+i] = 2.0*x[N] - x[N-i]
    end

    Lx = length(x_star)
    y = Vector{Float64}(undef, Lx - 1)
    @inbounds @simd for i in 1:Lx-1
        y[i] = (x_star[i+1] - x_star[i]) / tau0
    end
    Ly = length(y)
    cs = pushfirst!(cumsum(y), 0.0)

    for (k, m) in enumerate(m_values)
        nsubs = Ly - 3m + 1
        if nsubs < 1
            devs[k] = NaN
            continue
        end
        outer_sum = 0.0
        @inbounds @simd for j in 0:(nsubs - 1)
            h1 = cs[j+m+1]   - cs[j+1]
            h2 = cs[j+2m+1]  - cs[j+m+1]
            h3 = cs[j+3m+1]  - cs[j+2m+1]
            d3 = (h3 - 2.0*h2 + h1) / m
            outer_sum += d3^2
        end
        devs[k] = sqrt(outer_sum / (6.0 * nsubs))
    end

    return devs
end
```

Update dispatcher:

```julia
function _htotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:greenhall)
    if detrend === :greenhall || detrend === :legacy
        return _htotdev_greenhall(x, m_values, tau0)
    end
    detrend === :howe && return _htotdev_howe(x, m_values, tau0)
    throw(ArgumentError("unknown detrend recipe: $detrend; valid: :howe, :greenhall, :linear, :legacy"))
end
```

- [ ] **Step 4: Run, expect pass**

- [ ] **Step 5: Commit**

```
git commit -m "feat: HTOTDEV :howe recipe"
```

### Task 3.4: HTOTDEV `:linear`

**Files:**
- Modify: `lib/SigmaTauStability/src/core/total.jl`
- Modify: `lib/SigmaTauStability/test/runtests.jl`

- [ ] **Step 1: Smoke test**

```julia
@testset "HTOTDEV :linear smoke" begin
    using Random
    Random.seed!(20260508)
    N    = 1024
    tau0 = 1.0
    ms   = [1, 2, 4, 8, 16]
    x = _gen_powerlaw_phase(0.0, N; tau0=tau0)

    devs = SigmaTauStability._htotdev_core(x, ms, tau0; detrend=:linear)
    @test length(devs) == length(ms)
    @test all(isfinite, devs)
    @test all(>(0), devs)
    devs_legacy = SigmaTauStability._htotdev_core(x, ms, tau0; detrend=:legacy)
    @test all(0.1 .<= devs ./ devs_legacy .<= 10.0)
end
```

- [ ] **Step 2: Run, expect fail**

- [ ] **Step 3: Implement `_htotdev_linear`**

Derive from `_htotdev_greenhall` by replacing the half-mean slope estimate with full LS on the *frequency* series. Same time-reverse extension + third-diff operator.

```julia
function _htotdev_linear(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    # Per-window full LS detrend on the frequency series + time-reverse extension
    # + third-difference operator.
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    y = Vector{Float64}(undef, N-1)
    @inbounds @simd for i in 1:N-1
        y[i] = (x[i+1] - x[i]) / tau0
    end
    Ny = length(y)

    max_m = isempty(m_values) ? 0 : maximum(m_values)
    max_seg = 3 * max_m
    ext = Vector{Float64}(undef, 3 * max_seg)
    cs = Vector{Float64}(undef, 3 * max_seg + 1)

    for (k, m) in enumerate(m_values)
        if m == 1
            L = N - 3
            if L <= 0
                devs[k] = NaN
                continue
            end
            sum_sq = 0.0
            @inbounds @simd for i in 1:L
                d3 = x[i+3] - 3.0*x[i+2] + 3.0*x[i+1] - x[i]
                sum_sq += d3^2
            end
            devs[k] = sqrt(sum_sq / (6.0 * L * tau0^2))
            continue
        end

        n_iter = Ny - 3m + 1
        if n_iter < 1
            devs[k] = NaN
            continue
        end

        seg_len = 3m
        L_float = Float64(seg_len)
        sum_i = (L_float * (L_float + 1.0)) / 2.0
        sum_i2 = (L_float * (L_float + 1.0) * (2.0*L_float + 1.0)) / 6.0
        delta = L_float * sum_i2 - sum_i^2

        dev_sum = 0.0

        for i in 0:(n_iter - 1)
            sum_y = 0.0
            sum_iy = 0.0
            @inbounds @simd for j in 1:seg_len
                v = y[i+j]
                sum_y += v
                sum_iy += j * v
            end

            a = (sum_y * sum_i2 - sum_iy * sum_i) / delta
            b = (L_float * sum_iy - sum_y * sum_i) / delta

            @inbounds for j in 1:seg_len
                val = y[i+j] - (a + b * j)
                rev_val = y[i + seg_len - j + 1] - (a + b * (seg_len - j + 1))

                ext[j] = rev_val
                ext[seg_len + j] = val
                ext[2seg_len + j] = rev_val
            end

            cs[1] = 0.0
            @inbounds for j in 1:3seg_len
                cs[j+1] = cs[j] + ext[j]
            end

            sq = 0.0
            @inbounds @simd for j in 0:(6m - 1)
                h1 = (cs[j+m+1]  - cs[j+1])
                h2 = (cs[j+2m+1] - cs[j+m+1])
                h3 = (cs[j+3m+1] - cs[j+2m+1])
                sq += ((h3 - 2.0*h2 + h1) / m)^2
            end
            dev_sum += sq / (6.0 * m)
        end

        devs[k] = sqrt(dev_sum / (6.0 * n_iter))
    end

    return devs
end
```

Update dispatcher:

```julia
function _htotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:greenhall)
    if detrend === :greenhall || detrend === :legacy
        return _htotdev_greenhall(x, m_values, tau0)
    end
    detrend === :howe && return _htotdev_howe(x, m_values, tau0)
    detrend === :linear && return _htotdev_linear(x, m_values, tau0)
    throw(ArgumentError("unknown detrend recipe: $detrend; valid: :howe, :greenhall, :linear, :legacy"))
end
```

- [ ] **Step 4: Run, expect pass**

- [ ] **Step 5: Commit**

```
git commit -m "feat: HTOTDEV :linear recipe"
```

### Task 3.5: MHTOTDEV `:howe`

**Files:**
- Modify: `lib/SigmaTauStability/src/core/total.jl`
- Modify: `lib/SigmaTauStability/test/runtests.jl`

- [ ] **Step 1: Smoke test**

```julia
@testset "MHTOTDEV :howe smoke" begin
    using Random
    Random.seed!(20260508)
    N    = 1024
    tau0 = 1.0
    ms   = [1, 2, 4, 8]
    x = _gen_powerlaw_phase(0.0, N; tau0=tau0)

    devs = SigmaTauStability._mhtotdev_core(x, ms, tau0; detrend=:howe)
    @test length(devs) == length(ms)
    @test all(isfinite, devs)
    @test all(>(0), devs)
    devs_legacy = SigmaTauStability._mhtotdev_core(x, ms, tau0; detrend=:legacy)
    @test all(0.1 .<= devs ./ devs_legacy .<= 10.0)
end
```

- [ ] **Step 2: Run, expect fail**

- [ ] **Step 3: Implement `_mhtotdev_howe`**

MHTOT operator is averaged third-difference on phase. `:howe` does global mean-flip reflection on phase (no detrend), then applies the averaged third-diff.

```julia
function _mhtotdev_howe(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    # No detrend, global mean-flip reflection on phase. Apply averaged
    # third-difference operator on the extended series.
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    x_star = Vector{Float64}(undef, 3N - 4)
    @inbounds for i in 1:N-2
        x_star[i] = 2.0*x[1] - x[i+1]
    end
    @inbounds for i in 1:N
        x_star[N-2+i] = x[i]
    end
    @inbounds for i in 1:N-2
        x_star[2N-2+i] = 2.0*x[N] - x[N-i]
    end

    L = length(x_star)
    cs = pushfirst!(cumsum(x_star), 0.0)

    for (k, m) in enumerate(m_values)
        # Averaged third-difference: d3_j = ext[j] - 3 ext[j+m] + 3 ext[j+2m] - ext[j+3m]
        # then average over m successive starts (S[j+m] - S[j]) / m form via prefix sum.
        L3 = L - 3m
        if L3 < 1
            devs[k] = NaN
            continue
        end

        d3_vec = Vector{Float64}(undef, L3)
        @inbounds for j in 1:L3
            d3_vec[j] = x_star[j] - 3.0*x_star[j+m] + 3.0*x_star[j+2m] - x_star[j+3m]
        end

        S = pushfirst!(cumsum(d3_vec), 0.0)
        n_avg = L3 + 1 - m
        if n_avg <= 0
            devs[k] = NaN
            continue
        end

        block_var = 0.0
        @inbounds @simd for j in 1:n_avg
            block_var += (S[j+m] - S[j])^2
        end
        block_var /= (n_avg * 6.0 * Float64(m)^2)

        devs[k] = sqrt(block_var / (Float64(m)^2 * tau0^2))
    end

    return devs
end
```

Update dispatcher:

```julia
function _mhtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:linear)
    if detrend === :linear || detrend === :legacy
        return _mhtotdev_linear(x, m_values, tau0)
    end
    detrend === :howe && return _mhtotdev_howe(x, m_values, tau0)
    throw(ArgumentError("unknown detrend recipe: $detrend; valid: :howe, :greenhall, :linear, :legacy"))
end
```

- [ ] **Step 4: Run, expect pass**

- [ ] **Step 5: Commit**

```
git commit -m "feat: MHTOTDEV :howe recipe"
```

### Task 3.6: MHTOTDEV `:greenhall` + 5-noise-type smoke

**Files:**
- Modify: `lib/SigmaTauStability/src/core/total.jl`
- Modify: `lib/SigmaTauStability/test/runtests.jl`

- [ ] **Step 1: Smoke test (mid-spectrum + 5-noise-type)**

```julia
@testset "MHTOTDEV :greenhall smoke" begin
    using Random
    N    = 1024
    tau0 = 1.0
    ms   = [1, 2, 4, 8]

    # WHFM mid-spectrum check
    Random.seed!(20260508)
    x = _gen_powerlaw_phase(0.0, N; tau0=tau0)
    devs = SigmaTauStability._mhtotdev_core(x, ms, tau0; detrend=:greenhall)
    @test length(devs) == length(ms)
    @test all(isfinite, devs)
    @test all(>(0), devs)

    # 5-noise-type finite-output smoke (the new default needs basic coverage)
    for alpha in (2.0, 1.0, 0.0, -1.0, -2.0)
        Random.seed!(20260508)
        xa = _gen_powerlaw_phase(alpha, N; tau0=tau0)
        d = SigmaTauStability._mhtotdev_core(xa, ms, tau0; detrend=:greenhall)
        @test length(d) == length(ms)
        @test all(isfinite, d)
        @test all(>(0), d)
    end
end
```

- [ ] **Step 2: Run, expect fail**

- [ ] **Step 3: Implement `_mhtotdev_greenhall`**

Derive from `_mhtotdev_linear` by replacing the per-window full LS slope-and-intercept fit with a half-mean slope estimate. Same time-reverse extension + averaged third-diff.

```julia
function _mhtotdev_greenhall(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    # Per-window half-mean slope removal + time-reverse extension + averaged
    # third-difference operator. Same overall structure as _mhtotdev_linear,
    # but slope estimated by half-mean instead of full LS.
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    max_m = isempty(m_values) ? 0 : maximum(m_values)
    max_Lp = 3 * max_m + 1
    ext_len = 3 * max_Lp
    L3_max = ext_len - 3 * max_m
    ext = Vector{Float64}(undef, ext_len)
    d3_vec = Vector{Float64}(undef, L3_max)
    S = Vector{Float64}(undef, L3_max + 1)

    for (k, m) in enumerate(m_values)
        if m < 1
            devs[k] = NaN
            continue
        end
        nsubs = N - 4m + 1
        if nsubs < 1
            devs[k] = NaN
            continue
        end

        Lp = 3m + 1
        L3 = 3Lp - 3m

        total_sum = 0.0
        for n in 1:nsubs
            # Half-mean slope estimate over the Lp-length window
            half = floor(Int, Lp / 2)
            s1 = 0.0
            @inbounds @simd for j in 1:half
                s1 += x[n-1+j]
            end
            s1 /= half

            s2 = 0.0
            @inbounds @simd for j in (half+1):Lp
                s2 += x[n-1+j]
            end
            s2 /= (Lp - half)

            slope = (s2 - s1) / ((Lp / 2.0) * tau0)

            @inbounds for j in 1:Lp
                val = x[n-1+j] - slope * tau0 * (j - 1)
                rev_val = x[n-1 + Lp - j + 1] - slope * tau0 * (Lp - j)

                ext[j] = rev_val
                ext[Lp + j] = val
                ext[2Lp + j] = rev_val
            end

            @inbounds for j in 1:L3
                d3_vec[j] = ext[j] - 3.0*ext[j+m] + 3.0*ext[j+2m] - ext[j+3m]
            end

            S[1] = 0.0
            @inbounds for j in 1:L3
                S[j+1] = S[j] + d3_vec[j]
            end

            n_avg = L3 + 1 - m
            if n_avg > 0
                block_var = 0.0
                @inbounds @simd for j in 1:n_avg
                    block_var += (S[j+m] - S[j])^2
                end
                block_var /= (n_avg * 6.0 * Float64(m)^2)
            else
                block_var = 0.0
            end

            total_sum += block_var
        end

        devs[k] = sqrt(total_sum / (nsubs * Float64(m)^2 * tau0^2))
    end

    return devs
end
```

Update dispatcher:

```julia
function _mhtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:linear)
    if detrend === :linear || detrend === :legacy
        return _mhtotdev_linear(x, m_values, tau0)
    end
    detrend === :howe && return _mhtotdev_howe(x, m_values, tau0)
    detrend === :greenhall && return _mhtotdev_greenhall(x, m_values, tau0)
    throw(ArgumentError("unknown detrend recipe: $detrend; valid: :howe, :greenhall, :linear, :legacy"))
end
```

- [ ] **Step 4: Run, expect pass**

- [ ] **Step 5: Commit**

```
git commit -m "feat: MHTOTDEV :greenhall recipe + 5-noise-type smoke coverage"
```

**Phase 3 complete.** Push:

```
git push origin track-a4-detrend
```

This is a natural PR-B boundary if splitting (PR-A + PR-B = all infrastructure + new recipes, no defaults changed).

---

## Phase 4 — Default switch + breaking change polish

Switch TOTDEV default from `:legacy` to `:howe` and MHTOTDEV from `:linear` to `:greenhall`. Update parity tests to opt into `:legacy` explicitly. Tighten Stable32 TOTDEV cross-val. Add cross-recipe equivalence test. Update CHANGELOG and TODO.md.

### Task 4.1: Add cross-recipe equivalence test

**Files:**
- Modify: `lib/SigmaTauStability/test/runtests.jl`

- [ ] **Step 1: Add the testset**

```julia
@testset "Cross-recipe equivalence (legacy aliases)" begin
    # Verifies the alias claim in the design spec: for MTOT, HTOT, MHTOT,
    # :legacy must be a strict alias for :greenhall (MTOT/HTOT) or :linear
    # (MHTOT). Tested at rtol=1e-12 on the legacy_kernels.jl fixture.
    using Random
    Random.seed!(20260507)
    N    = 4096
    tau0 = 1.0
    wpm  = randn(N) .* 1e-9
    rwfm = cumsum(cumsum(randn(N) .* 1e-12))
    x    = wpm .+ rwfm

    rt = 1e-12
    at = 1e-25

    for m in [1, 2, 4, 8, 16]
        @test SigmaTauStability._mtotdev_core(x, [m], tau0; detrend=:legacy)[1]    ≈
              SigmaTauStability._mtotdev_core(x, [m], tau0; detrend=:greenhall)[1] atol=at rtol=rt
        @test SigmaTauStability._htotdev_core(x, [m], tau0; detrend=:legacy)[1]    ≈
              SigmaTauStability._htotdev_core(x, [m], tau0; detrend=:greenhall)[1] atol=at rtol=rt
    end

    for m in [1, 2, 4, 8]
        @test SigmaTauStability._mhtotdev_core(x, [m], tau0; detrend=:legacy)[1] ≈
              SigmaTauStability._mhtotdev_core(x, [m], tau0; detrend=:linear)[1] atol=at rtol=rt
    end
end
```

- [ ] **Step 2: Run, expect pass**

```
pkg> test SigmaTauStability
```

- [ ] **Step 3: Commit**

```
git add lib/SigmaTauStability/test/runtests.jl
git commit -m "test: cross-recipe equivalence (legacy aliases for MTOT/HTOT/MHTOT)"
```

### Task 4.2: Switch TOTDEV default to `:howe`

**Files:**
- Modify: `lib/SigmaTauStability/src/core/total.jl` (one-line default change)
- Modify: `lib/SigmaTauStability/src/api/total.jl` (default for `totdev`)
- Modify: `lib/SigmaTauStability/test/runtests.jl` (legacy parity TOTDEV calls + Stable32 floor)

- [ ] **Step 1: Update legacy_kernels.jl parity test calls**

In `runtests.jl`, find the `Legacy parity (extracted SP1065 kernels)` testset (around line 114). Update the TOTDEV loop to pass `detrend=:legacy` explicitly:

```julia
# TOTDEV. (Smaller grid — kernel is O(N) per m but allocates an extended
# 3N-4 array each call.)
for m in [1, 2, 4, 8, 16, 32]
    new_dev = sqrt(LK.totdev_var(x, m, tau0))
    @test SigmaTauStability._totdev_core(x, [m], tau0; detrend=:legacy)[1] ≈ new_dev atol=at rtol=rt
end
```

The other TOT-family parity tests (MTOT/HTOT/MHTOT) don't need `detrend=:legacy` because:
- MTOT/HTOT have `:greenhall` as default; current legacy_kernels.jl behavior IS `:greenhall`. They pass either way.
- MHTOT will be switched to `:greenhall` in Task 4.3; then *its* parity test will need `detrend=:legacy`.

- [ ] **Step 2: Switch TOTDEV default in core**

In `core/total.jl`:

```julia
function _totdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:howe)
    detrend === :howe && return _totdev_howe(x, m_values, tau0)
    detrend === :legacy && return _totdev_legacy(x, m_values, tau0)
    detrend === :greenhall && return _totdev_greenhall(x, m_values, tau0)
    detrend === :linear && return _totdev_linear(x, m_values, tau0)
    throw(ArgumentError("unknown detrend recipe: $detrend; valid: :howe, :greenhall, :linear, :legacy"))
end
```

- [ ] **Step 3: Switch TOTDEV default in API wrapper**

In `api/total.jl`:

```julia
function totdev(data::PhaseData, m_values::Vector{Int}; detrend::Symbol=:howe, calc_ci::Bool=true, confidence::Float64=0.95)
    ...
end
```

- [ ] **Step 4: Tighten Stable32 TOTDEV cross-val**

In `runtests.jl` find the `@testset "Stable32 cross-validation (reference/validation/)"` block (around line 207). Update the `"Total"` branch (around line 266) from rtol=0.15 to rtol=1e-4. The kernel call now uses `:howe` by default — the assertion can stay using `LK.totdev_var(...)` since legacy_kernels.jl is verbatim MATLAB-era reference and was kept for parity, OR we can switch to using the new `:howe` kernel directly.

Cleanest: update the `"Total"` branch to use the new core function with explicit `:howe`:

```julia
elseif kind == "Total"
    # SigmaTau :howe matches Stable32's TOTDEV (SP1065 eqn 25) at rtol=1e-4.
    # The :legacy global-LS detrend recipe diverges by ~15% at long τ; not
    # exercised here.
    got = SigmaTauStability._totdev_core(x, [m], tau0; detrend=:howe)[1]
    @test got ≈ sigma_ref rtol=1e-4
```

- [ ] **Step 5: Run, expect pass**

```
pkg> test SigmaTauStability
```

- [ ] **Step 6: Commit**

```
git add lib/SigmaTauStability/src/core/total.jl lib/SigmaTauStability/src/api/total.jl lib/SigmaTauStability/test/runtests.jl
git commit -m "feat!: TOTDEV default switches to :howe (SP1065 eqn 25)"
```

The `feat!` exclamation flags the breaking change per Conventional Commits.

### Task 4.3: Switch MHTOTDEV default to `:greenhall`

**Files:**
- Modify: `lib/SigmaTauStability/src/core/total.jl`
- Modify: `lib/SigmaTauStability/src/api/total.jl`
- Modify: `lib/SigmaTauStability/test/runtests.jl`

- [ ] **Step 1: Update legacy_kernels.jl parity for MHTOT**

In `runtests.jl`'s `Legacy parity` testset, update the MHTOTDEV loop:

```julia
# MHTOTDEV.
for m in [1, 2, 4, 8]
    new_dev = sqrt(LK.mhtotdev_var(x, m, tau0))
    @test SigmaTauStability._mhtotdev_core(x, [m], tau0; detrend=:legacy)[1] ≈ new_dev atol=at rtol=rt
end
```

- [ ] **Step 2: Switch MHTOTDEV default in core**

```julia
function _mhtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:greenhall)
    detrend === :greenhall && return _mhtotdev_greenhall(x, m_values, tau0)
    detrend === :linear && return _mhtotdev_linear(x, m_values, tau0)
    detrend === :legacy && return _mhtotdev_linear(x, m_values, tau0)  # alias
    detrend === :howe && return _mhtotdev_howe(x, m_values, tau0)
    throw(ArgumentError("unknown detrend recipe: $detrend; valid: :howe, :greenhall, :linear, :legacy"))
end
```

- [ ] **Step 3: Switch MHTOTDEV default in API wrapper**

```julia
function mhtotdev(data::PhaseData, m_values::Vector{Int}; detrend::Symbol=:greenhall, calc_ci::Bool=true, confidence::Float64=0.95)
    ...
end
```

- [ ] **Step 4: Run, expect pass**

```
pkg> test SigmaTauStability
```

- [ ] **Step 5: Commit**

```
git add lib/SigmaTauStability/src/core/total.jl lib/SigmaTauStability/src/api/total.jl lib/SigmaTauStability/test/runtests.jl
git commit -m "feat!: MHTOTDEV default switches to :greenhall (was :linear)"
```

### Task 4.4: CHANGELOG entry

**Files:**
- Modify: `CHANGELOG.md`

- [ ] **Step 1: Add Changed entry under `[Unreleased]`**

```markdown
- **Breaking:** `totdev` default detrend recipe is now `:howe` (SP1065 eqn 25:
  no detrend, mean-flip endpoint reflection). Previous behavior (global LS
  detrend on top of the same reflection) is available via
  `totdev(...; detrend=:legacy)`. Output values change for all τ. The new
  default matches Stable32 / allantools at `rtol=1e-4` (vs the prior
  `rtol=0.15` boundary-policy floor).
- **Breaking:** `mhtotdev` default detrend recipe is now `:greenhall`
  (per-window half-mean slope removal + time-reverse extension), matching
  the convention of MTOT / HTOT in the Hadamard-modified family. Previous
  behavior (per-window full LS detrend) is available via
  `mhtotdev(...; detrend=:legacy)` or equivalently `:linear`. MHTOTDEV is
  novel to SigmaTau; no external numerical reference exists, so the
  recipe choice is a methodology decision rather than a parity contract.
- New `detrend::Symbol` kwarg on all four total-family kernels and API
  wrappers (`totdev`, `mtotdev`, `htotdev`, `mhtotdev`). Recipes:
  `:howe`, `:greenhall`, `:linear`, `:legacy`. See spec at
  `docs/superpowers/specs/2026-05-07-detrend-kwarg-design.md` for the
  per-recipe math.
- `totdev`, `mtotdev`, `htotdev`, `mhtotdev` API surface tightened by
  removing the silent global-LS detrend in `totdev`'s default path. Bias
  correction `bias_correction(:totvar, ...)` is now correctly calibrated
  against the new `:howe` default (SP1065 formula).
- MTOTDEV, HTOTDEV default outputs are unchanged (their previous
  behavior was already `:greenhall`).
```

- [ ] **Step 2: Commit**

```
git add CHANGELOG.md
git commit -m "docs: CHANGELOG entry for detrend kwarg + default switches"
```

### Task 4.5: TODO.md follow-up entries

**Files:**
- Modify: `TODO.md`

- [ ] **Step 1: Add follow-up items**

In `TODO.md`, under the appropriate priority section (Medium for the Monte Carlo, Low for the cross-val tightening):

```markdown
- [ ] **MHTOTDEV bias / EDF Monte Carlo.** Synthesize known-noise via
  `_gen_powerlaw_phase` for each α; compute MHTOT and MHDEV; the ratio
  yields the bias factor B(α). Re-fit `_coeff_mhtot` empirically per
  detrend recipe (`:greenhall`, `:linear`). Track per-recipe — EDF is
  recipe-specific. Spec: `docs/superpowers/specs/2026-05-07-detrend-kwarg-design.md`
  → "Out-of-scope / future work".
- [ ] **allantools cross-validation TOTDEV tightening.** After
  `track-b1-allantools` merges to main, update
  `lib/SigmaTauStability/test/allantools_cross_validation.jl` so the
  `"Total"` branch uses `_totdev_core(...; detrend=:howe)` and rtol=1e-4
  (currently rtol=0.15 with `:legacy`). One-line change in that test
  file.
```

- [ ] **Step 2: Commit**

```
git add TODO.md
git commit -m "docs: TODO entries for MHTOT EDF Monte Carlo + allantools TOT tightening"
```

**Phase 4 complete.** Push:

```
git push origin track-a4-detrend
```

This is a natural PR-C boundary. The PR includes the breaking changes; CHANGELOG flags them clearly.

---

## Final verification

- [ ] **Step 1: Full test run**

```
pkg> test SigmaTauStability
```

Expected counts:
- Existing 339 assertions still pass (legacy parity now opts into `:legacy`)
- 7 new testsets: TOTDEV {`:howe` tight, `:greenhall` smoke, `:linear` smoke}, MTOTDEV {`:howe` smoke, `:linear` smoke}, HTOTDEV {`:howe` smoke, `:linear` smoke}, MHTOTDEV {`:howe` smoke, `:greenhall` smoke + 5-noise-type}, Cross-recipe equivalence
- ~25–30 new assertions

Total roughly 365–370 passing.

- [ ] **Step 2: Sanity-check the `using SigmaTau` smoke**

```
julia> using SigmaTau
julia> using Random; Random.seed!(7)
julia> data = PhaseData(randn(1000) .* 1e-9, 1.0)
julia> totdev(data, [1, 2, 4]; detrend=:howe)     # canonical SP1065 default
julia> totdev(data, [1, 2, 4]; detrend=:legacy)   # opt-in old behavior
julia> mhtotdev(data, [1, 2, 4]; detrend=:greenhall)  # new default
julia> mhtotdev(data, [1, 2, 4]; detrend=:linear)     # opt-in old behavior
```

All four calls should return finite, positive `StabilityResult` outputs.

- [ ] **Step 3: Open PR(s)**

If shipping as a single PR:

```
gh pr create --base main --head track-a4-detrend \
  --title "feat!: detrend kwarg for total-family kernels (Track A4)" \
  --body-file <(cat <<EOF
Adds a \`detrend::Symbol\` kwarg to totdev / mtotdev / htotdev / mhtotdev
with four named recipes (:howe, :greenhall, :linear, :legacy) and
canonical per-kernel defaults.

Breaking changes:
- totdev default: :legacy → :howe (SP1065 eqn 25, matches Stable32 at rtol=1e-4)
- mhtotdev default: :legacy(=:linear) → :greenhall

Spec: docs/superpowers/specs/2026-05-07-detrend-kwarg-design.md

## Test plan
- [x] All 339 existing tests pass with :legacy opt-in for TOT and MHTOT parity
- [x] New: TOTDEV :howe tight match against Stable32 at rtol=1e-4
- [x] New: per-recipe smokes for all 16 (recipe, kernel) cells
- [x] New: cross-recipe equivalence (MTOT/HTOT :legacy ≡ :greenhall)
EOF
)
```

If splitting into PR-A / PR-B / PR-C: open them sequentially against main, each one based on the previous, with the breaking-change PR clearly flagged.

---

## Self-review notes

- All 4 spec recipes covered: ✓ `:howe` (Tasks 2.1, 3.1, 3.3, 3.5), `:greenhall` (Task 2.2, also via aliases for MTOT/HTOT, and Task 3.6 for MHTOT), `:linear` (Tasks 2.3, 3.2, 3.4, also via alias for MHTOT), `:legacy` (Phase 1 refactors).
- All 4 spec defaults covered: ✓ TOTDEV → `:howe` (Task 4.2), MTOT/HTOT → `:greenhall` (Phase 1, no behavior change), MHTOT → `:greenhall` (Task 4.3).
- Spec Section 4 testing strategy covered:
  - Legacy parity preserved: ✓ Tasks 4.2 step 1 and 4.3 step 1 update the two parity calls that need explicit `:legacy`.
  - Stable32 TOTDEV tightening: ✓ Task 4.2 step 4.
  - Per-recipe smokes for all 16 cells: ✓ Tasks 2.1–2.3, 3.1–3.6 cover the 9 new cells; the 7 existing-as-default cells are covered by the existing test corpus.
  - MHTOT 5-noise-type smoke: ✓ Task 3.6 step 1.
  - Cross-recipe equivalence (`:legacy` aliases): ✓ Task 4.1.
  - allantools tightening: tracked as TODO in Task 4.5 (deferred per spec).
- No placeholders.
- Type / function name consistency: dispatcher signatures and helper names match across phases.
