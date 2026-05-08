# SigmaTauEnsemble API

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

## Estimators

```@docs
AbstractEstimator
StandardKalmanFilter
UDFactorizedFilter
KuramotoOscillator
predict!
update!
```

## Steering

```@docs
PIDController
step!
steer_to_correction
```
