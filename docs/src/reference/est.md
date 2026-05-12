# SigmaTau.Est — Estimation

## Clock models

```@docs
AbstractClockModel
TwoStateClock
ThreeStateClock
RelativisticClock
nstates
state_transition
process_noise
measurement_matrix
measurement_noise
```

## Clock ensembles

```@docs
ClockEnsemble
EnsembleWeights
```

## Estimators

```@docs
AbstractEstimator
StandardKalmanFilter
UDFactorizedFilter
KuramotoOscillator
predict!
update!
prop!
```

## Steering

```@docs
PIDController
step!
steer_to_correction
```
