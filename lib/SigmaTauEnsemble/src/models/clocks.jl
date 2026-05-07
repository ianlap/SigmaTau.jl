# models/clocks.jl — Clock Model Definitions

abstract type AbstractClockModel end

Base.@kwdef struct TwoStateClock <: AbstractClockModel
    tau::Float64
    q0::Float64 = 0.0 # WPM (measurement)
    q1::Float64 = 0.0 # WFM (state)
    q2::Float64 = 0.0 # RWFM (state)
end

Base.@kwdef struct ThreeStateClock <: AbstractClockModel
    tau::Float64
    q0::Float64 = 0.0 # WPM (measurement)
    q1::Float64 = 0.0 # WFM (state)
    q2::Float64 = 0.0 # RWFM (state)
    q3::Float64 = 0.0 # IRWFM (state)
end

struct RelativisticClock <: AbstractClockModel end # Lunar PNT specifics

nstates(::TwoStateClock) = 2
nstates(::ThreeStateClock) = 3

function state_transition(m::TwoStateClock)
    @SMatrix [1.0 m.tau; 0.0 1.0]
end

function state_transition(m::ThreeStateClock)
    @SMatrix [1.0 m.tau m.tau^2 / 2.0; 0.0 1.0 m.tau; 0.0 0.0 1.0]
end

function process_noise(m::TwoStateClock)
    τ = m.tau
    Q11 = m.q1*τ + m.q2*τ^3/3.0
    Q12 = m.q2*τ^2/2.0
    Q22 = m.q2*τ
    @SMatrix [Q11 Q12; Q12 Q22]
end

function process_noise(m::ThreeStateClock)
    τ = m.tau
    τ2 = τ^2; τ3 = τ^3; τ4 = τ^4; τ5 = τ^5
    Q11 = m.q1*τ + m.q2*τ3/3.0 + m.q3*τ5/20.0
    Q12 = m.q2*τ2/2.0 + m.q3*τ4/8.0
    Q13 = m.q3*τ3/6.0
    Q22 = m.q2*τ + m.q3*τ3/3.0
    Q23 = m.q3*τ2/2.0
    Q33 = m.q3*τ
    @SMatrix [Q11 Q12 Q13; Q12 Q22 Q23; Q13 Q23 Q33]
end

measurement_matrix(::TwoStateClock) = @SMatrix [1.0 0.0]
measurement_matrix(::ThreeStateClock) = @SMatrix [1.0 0.0 0.0]

measurement_noise(m::TwoStateClock) = @SMatrix [m.q0]
measurement_noise(m::ThreeStateClock) = @SMatrix [m.q0]
