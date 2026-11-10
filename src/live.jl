# Live measurement: a source gives one reading at a time, and a set of
# streaming accumulators turns the readings into running deviations.
#
# A source is anything with `read_reading(src) → Float64`, NaN once the
# stream has ended. `InstrumentSource` reads lines of text from an IO: a
# serial device opened as a file, a TCP socket to an instrument (SCPI on
# port 5025) or a pipe. `SimulatedSource` is a noise model in its place.

"""
    read_reading(src) → Float64

Next reading of a live source, or `NaN` once its stream has ended. A source
is a [`SimulatedSource`](@ref), an [`InstrumentSource`](@ref), or any type
given a method of this function.
"""
function read_reading end

"""
    SimulatedSource(; white_pm=0, white_fm=1e-12, rw_fm=0, tau0=1.0, rng=nothing)

Simulated clock, for trying a live setup without an instrument. A reading
is a phase sample in seconds, the sum of three independent Gaussian terms.

# Keywords

- `white_pm=0`: standard deviation of the white PM added to each phase
  sample, in seconds.
- `white_fm=1e-12`: standard deviation of the white FM, as a fractional
  frequency per sample.
- `rw_fm=0`: standard deviation of the step the frequency takes from one
  sample to the next, which makes random-walk FM.
- `tau0=1.0`: sample interval in seconds.
- `rng`: an `AbstractRNG`, for a reproducible stream.

At ``τ_0`` the Allan deviation of the white PM term is
`√3·white_pm/tau0`, that of the white FM term `white_fm`, and that of the
random-walk FM term `rw_fm/√2`.

# Examples

```julia
clock = SimulatedSource(white_fm=1e-12)
live = LiveStability([:adev], 1.0, [1, 10, 100])
for _ in 1:10_000
    push!(live, read_reading(clock))
end
snapshot(live)[1]
```
"""
mutable struct SimulatedSource{R}
    const white_pm::Float64
    const white_fm::Float64
    const rw_fm::Float64
    const tau0::Float64
    const rng::R
    x::Float64          # integrated frequency, seconds
    y_rw::Float64       # random-walk frequency
end

function SimulatedSource(; white_pm::Real=0.0, white_fm::Real=1e-12, rw_fm::Real=0.0,
                         tau0::Real=1.0, rng=nothing)
    tau0 > 0 || throw(ArgumentError("SimulatedSource: tau0 must be positive"))
    min(white_pm, white_fm, rw_fm) >= 0 || throw(ArgumentError(
        "SimulatedSource: noise levels must be non-negative"))
    return SimulatedSource(Float64(white_pm), Float64(white_fm), Float64(rw_fm),
                           Float64(tau0), rng, 0.0, 0.0)
end

_draw(::Nothing) = randn()
_draw(rng) = randn(rng)

function read_reading(s::SimulatedSource)
    pm = s.white_pm * _draw(s.rng)
    s.y_rw += s.rw_fm * _draw(s.rng)
    s.x += (s.white_fm * _draw(s.rng) + s.y_rw) * s.tau0
    return s.x + pm
end

"""
    InstrumentSource(io; query="")

Readings from an instrument that sends lines of text. A reading is the
first number on a line, as in a SCPI reply such as `+1.000000000012E+07`
or in comma-separated fields. Blank lines are skipped, and the stream ends
with `io`.

With a `query` such as `"READ?"`, the query is written before each line is
read. Without one the instrument is expected to send readings on its own.

SCPI's not-a-number value, `9.91E+37`, with which a counter reports a gate
without a measurement, throws an error: the running sums cannot leave a
sample out.

`io` is any open stream:

- an instrument on the network, `Sockets.connect(host, 5025)`: 5025 is the
  port registered for raw SCPI;
- a serial port, opened as a file once its line speed is set: on Linux
  `stty -F /dev/ttyUSB0 9600 raw`, on macOS `stty -f`, then
  `open("/dev/ttyUSB0", "r+")`;
- a pipe, or a file of readings, read to its end.

# Examples

```julia
using Sockets
counter = InstrumentSource(Sockets.connect("192.168.1.20", 5025); query="READ?")
live = LiveStability([:adev, :mdev], 1.0, [1, 2, 4, 8, 16];
                     quantity=:hz, nominal=10e6)
for _ in 1:3600
    push!(live, read_reading(counter))
end
snapshot(live)
```
"""
struct InstrumentSource{T<:IO}
    io::T
    query::String
end
InstrumentSource(io::IO; query::AbstractString="") = InstrumentSource(io, String(query))

function read_reading(s::InstrumentSource)
    if !isempty(s.query)
        write(s.io, s.query, '\n')
        flush(s.io)
    end
    while !eof(s.io)
        fields = split(readline(s.io), (',', ';', ' ', '\t', '\r'); keepempty=false)
        isempty(fields) && continue
        v = tryparse(Float64, fields[1])
        v === nothing && throw(ArgumentError("InstrumentSource: not a reading: \"" * fields[1] * "\""))
        # A missed gate is a gap, and a running estimate cannot bridge one.
        isnan(_scpi_nan(v)) && throw(ArgumentError(
            "InstrumentSource: the instrument returned no reading (9.91E+37, a missed gate)"))
        return v
    end
    return NaN
end

"""
    LiveStability(devs, tau0, m_values; quantity=:phase, nominal=NaN)

Running deviations of a live measurement: one [`StreamingStability`](@ref)
for each name in `devs`, all fed by `push!(live, reading)`.
[`snapshot`](@ref) returns one [`StabilityResult`](@ref) per deviation,
equal to the estimator's on the readings so far, and [`nsamples`](@ref) the
number of readings.

# Arguments

- `devs`: a vector of `:adev`, `:mdev`, `:hdev` and `:mhdev`.
- `tau0`: interval between readings, in seconds.
- `m_values`: averaging factors.

# Keywords

- `quantity=:phase`: what a reading is: `:phase`, a time error in seconds;
  `:frequency`, a fractional frequency, integrated to phase as
  [`FrequencyData`](@ref) is, so that `n` readings make `n + 1` phase
  samples; or `:hz`, a frequency in hertz as a counter reads it, converted
  to fractional frequency against `nominal`.
- `nominal`: the nominal frequency in hertz for `:hz`. By default it is the
  first reading rounded to four significant digits.

A reading that is not finite throws an error.
"""
mutable struct LiveStability
    const accs::Vector{StreamingStability}
    const quantity::Symbol
    const tau0::Float64
    nominal::Float64
    x::Float64
    readings::Int
end

function LiveStability(devs::AbstractVector{Symbol}, tau0::Real,
                       m_values::AbstractVector{<:Integer};
                       quantity::Symbol=:phase, nominal::Real=NaN)
    quantity in (:phase, :frequency, :hz) || throw(ArgumentError(
        "LiveStability: quantity must be :phase, :frequency or :hz"))
    isempty(devs) && throw(ArgumentError("LiveStability: devs must be non-empty"))
    accs = StreamingStability[StreamingStability(d, tau0, m_values) for d in devs]
    return LiveStability(accs, quantity, Float64(tau0), Float64(nominal), 0.0, 0)
end

function Base.push!(live::LiveStability, reading::Real)
    r = Float64(reading)
    isfinite(r) || throw(ArgumentError("LiveStability: readings must be finite"))
    q = live.quantity
    if q === :phase
        _push_phase!(live, r)
    else
        if q === :hz
            isnan(live.nominal) && (live.nominal = round(r; sigdigits=4))
            r = r / live.nominal - 1.0
        end
        live.readings == 0 && _push_phase!(live, 0.0)
        live.x += r * live.tau0
        _push_phase!(live, live.x)
    end
    live.readings += 1
    return live
end

function _push_phase!(live::LiveStability, x::Float64)
    for acc in live.accs
        push!(acc, x)
    end
end

snapshot(live::LiveStability) = StabilityResult[snapshot(acc) for acc in live.accs]

"""
    nsamples(live::LiveStability) → Int

Number of readings received.
"""
nsamples(live::LiveStability) = live.readings
