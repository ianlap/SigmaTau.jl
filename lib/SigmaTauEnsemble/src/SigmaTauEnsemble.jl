module SigmaTauEnsemble

using DocStringExtensions
using LinearAlgebra
using StaticArrays

include("models/clocks.jl")
include("estimators/filters.jl")

export AbstractClockModel, TwoStateClock, ThreeStateClock, RelativisticClock
export nstates, state_transition, process_noise, measurement_matrix, measurement_noise
export AbstractEstimator, StandardKalmanFilter, UDFactorizedFilter, KuramotoOscillator
export predict!, update!
export PIDController, step!, steer_to_correction

end # module
