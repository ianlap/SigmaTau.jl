# Power-law noise

The random part of a clock's instability is modeled as a sum of power laws
in the spectrum of its fractional frequency
([SP1065 §3.2](@cite riley-2008-sp1065);
[IEEE Std 1139-2022, Annex B](@cite ieee1139-2022-definitions)):

```math
S_y(f) = \sum_α h_α f^α .
```

Each term is a noise type with its own slope on a stability plot. That is
what lets a plot be read, and it is what every confidence interval and bias
factor in the package depends on.

## The noise types

| Label | Noise | ``α`` | ADEV, HDEV | MDEV, MHDEV | TDEV, HTDEV |
|:---|:---|---:|:---|:---|:---|
| `:WHPM` | white PM | 2 | ``τ^{-1}`` | ``τ^{-3/2}`` | ``τ^{-1/2}`` |
| `:FLPM` | flicker PM | 1 | ``τ^{-1}``, nearly | ``τ^{-1}`` | ``τ^{0}`` |
| `:WHFM` | white FM | 0 | ``τ^{-1/2}`` | ``τ^{-1/2}`` | ``τ^{1/2}`` |
| `:FLFM` | flicker FM | −1 | ``τ^{0}`` | ``τ^{0}`` | ``τ^{1}`` |
| `:RWFM` | random-walk FM | −2 | ``τ^{1/2}`` | ``τ^{1/2}`` | ``τ^{3/2}`` |
| `:FWFM` | flicker-walk FM | −3 | ``τ^{1}``, HDEV only | ``τ^{1}``, MHDEV only | ``τ^{2}``, HTDEV only |
| `:RRFM` | random-run FM | −4 | ``τ^{3/2}``, HDEV only | ``τ^{3/2}``, MHDEV only | ``τ^{5/2}``, HTDEV only |

The columns give the slope of each deviation on log-log axes. The two phase
noises share the slope ``τ^{-1}`` under the Allan deviation, flicker PM with
a logarithmic factor, and the modified deviations separate them
([SP1065 §5.1](@cite riley-2008-sp1065)). The Allan variances do not
converge for the last two rows; the third difference of the Hadamard family
does ([SP1065 §3.4](@cite riley-2008-sp1065)). The total and Thêo estimators
follow the slope of the deviation they extend.

A linear frequency drift ``D``, in fractional frequency per second, is not
noise. Its phase ``Dt^2/2`` has the second difference ``Dτ^2``, which adds
``D^2τ^2/2`` to the Allan variance: the deviation rises as ``τ^{1}`` at
long averaging times ([SP1065 §5.2](@cite riley-2008-sp1065)). The third
difference is zero, and the drift adds nothing to the Hadamard family.

```@example noise
using SigmaTau, Plots, Random

names = ["white PM", "flicker PM", "white FM", "flicker FM", "random-walk FM"]
plt = plot(layout = (1, 2), size = (760, 330), legend = :bottomleft)
for (k, α) in enumerate(2:-1:-2)
    p = noise_gen(PhaseData, 16384, 1.0; sigma1 = Dict(α => 1e-12), rng = Xoshiro(α + 3))
    plot!(plt, adev(p; ci = false); subplot = 1, label = names[k], title = "ADEV")
    plot!(plt, mdev(p; ci = false); subplot = 2, label = "", title = "MDEV")
end
plt
```

Each record has ``σ_y(τ_0) = 10^{-12}``. White and flicker PM run together
on the left and apart on the right.

## From a spectrum to a deviation

A term ``h_α f^α`` contributes a known variance to each estimator. For the
Allan variance ([SP1065 §7.1, Eq. (67)](@cite riley-2008-sp1065);
[IEEE Std 1139-2022, Table B.2](@cite ieee1139-2022-definitions)):

| Noise | ``σ_y^2(τ)`` |
|:---|:---|
| white PM | ``3 f_h h_2 / (4π^2 τ^2)`` |
| flicker PM | ``[1.038 + 3\ln(2π f_h τ)]\, h_1 / (4π^2 τ^2)`` |
| white FM | ``h_0 / (2τ)`` |
| flicker FM | ``2\ln 2 \; h_{-1}`` |
| random-walk FM | ``(2π^2/3)\, h_{-2}\, τ`` |

The two PM terms depend on the measurement bandwidth ``f_h``.

These forms describe a continuous signal observed through a sharp low-pass
filter. A record is a sampled signal, whose spectrum folds at the Nyquist
frequency ``1/(2τ_0)``, and the two differ at the shortest averaging times.
[`powerlaw_sigma`](@ref) integrates the spectrum of the sampled record
through the transfer function of the estimator, for ADEV, MDEV, HDEV and
MHDEV, and so is exact from ``τ_0``. On flicker FM the sampled record lies
44 % above the continuous form in variance at ``τ_0`` and on random-walk FM
50 %, both within 1 % from ``16τ_0``. On flicker PM it is 19 % above at
``τ_0`` and still 4 % above at ``64τ_0``. The white PM and white FM forms
hold at every ``τ``.

```@example noise
h = Dict(:WHFM => 2e-24, :RWFM => 1e-30)       # S_y(f) = 2e-24 + 1e-30 f⁻²
round.(powerlaw_sigma(h, [1.0, 100.0, 1e4]; tau0 = 1.0); sigdigits = 4)
```

[`powerlaw_psd`](@ref) gives the same model as ``S_y(f)``, ``S_x(f)`` or
``ℒ(f)``, and [`powerlaw_fit`](@ref) goes the other way, from a measured
deviation or spectrum to the ``h_α``.

!!! note "Hadamard conversions"
    Table B.2 of IEEE Std 1139-2022 also has a column for the Hadamard
    variance. Its entries are two thirds of what the variance of the
    standard's own Eq. (A.31) gives, for example ``h_0/(3τ)`` on white FM
    where the definition gives ``h_0/(2τ)``. The package follows the
    definition, as [SP1065 §5.2.8](@cite riley-2008-sp1065) does.

## Simulated noise

[`noise_gen`](@ref) draws each component as white Gaussian noise through the
fractional differencing filter ``(1 - z^{-1})^{α/2}`` of
[Kasdin and Walter (1992)](@cite kasdin-1992-discrete-simulation). A
component asked for by its coefficient ``h_α`` has the spectrum

```math
S_y(f) = h_α \left(\frac{\sin πfτ_0}{πτ_0}\right)^{α},
```

which is ``h_α f^α`` well below the Nyquist frequency. For ``α ≤ -1`` the
process is not stationary and the spectrum is that of its differences. The
test suite averages the Welch spectrum of 40 records for each ``α`` from 2
to −4 and holds its mean ratio to this formula, across the band, within
3 %.

A component can also be asked for by its Allan deviation at ``τ_0``, which
is how most examples in this documentation are built:

```@example noise
p = noise_gen(PhaseData, 4096, 1.0; sigma1 = Dict(:WHFM => 1e-12, :RWFM => 3e-14),
              rng = Xoshiro(1))
adev(p, [1, 16, 256])
```

## Fitting a model to a record

Two functions decompose a measured curve into power laws.
[`powerlaw_fit`](@ref) returns the calibrated ``h_α`` by nonnegative least
squares, from a deviation or from a spectrum, and
[`noise_id_nnls`](@ref) names the noise types present. [`noise_fit`](@ref)
builds the same model by hand, one noise type over one range of ``τ`` at a
time, for a record whose curve a blind fit reads badly.

```@example noise
f = powerlaw_fit(p; alphas = [0, -2])
```

The coefficients cannot be negative, so a noise type the record does not
contain comes out small and positive. Compare `f.fit` with `f.y` before
relying on a small coefficient.
