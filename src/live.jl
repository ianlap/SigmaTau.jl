# Live measurement: instrument sources that yield one reading at a time, and
# a set of streaming accumulators that turn readings into running σ(τ).
#
# A source is anything with `read_reading(src) → Float64` (NaN once the
# stream ends). `InstrumentSource` reads text lines from any IO: a serial
# device opened as a file, a TCP socket to an LXI instrument (SCPI on port
# 5025), or a pipe. `SimulatedSource` is a noise-only stand-in clock.

"""
    read_reading(src) → Float64

The next reading from a live source, or `NaN` once the stream has ended.
"""
function read_reading end

"""
    SimulatedSource(; white_pm=0, white_fm=1e-12, rw_fm=0, tau0=1.0, rng=nothing)

A simulated clock for trying out a live setup without an instrument. Each
reading is a phase sample in seconds, built by exact recurrences from
independent Gaussian draws: white PM adds `white_pm` seconds (standard
deviation) to each phase sample, white FM is a fractional frequency of
standard deviation `white_fm` per sample, and random-walk FM steps that
frequency by `rw_fm` per sample. At τ₀ the expected Allan deviations of
the three terms are `white_fm`, `√3·white_pm/tau0` and `rw_fm/√2`.

Pass an `AbstractRNG` as `rng` for a reproducible stream.
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
    InstrumentSource(io::IO; query="")

Readings from an instrument that talks in text lines: each reading is the
first number on a line (SCPI replies such as `+1.000000000012E+07`, or
comma-separated fields). With a `query` such as `"READ?"`, the query is
written before every line is read; without one the instrument is expected
to stream on its own. Blank lines are skipped and the stream ends at EOF.
SCPI's no-reading value `9.91E+37` (a missed gate) throws: a running
estimate cannot bridge a gap.

`io` is any open stream:

- an LXI instrument over the network: `Sockets.connect(host, 5025)` (raw
  SCPI socket; most counters and VISA "SOCKET" resources use port 5025);
- a serial port, opened as a file after setting its line speed, e.g.
  `stty -F /dev/ttyUSB0 9600 raw` (Linux; `stty -f` on macOS) and then
  `open("/dev/ttyUSB0", "r+")`;
- a pipe, or a log file, replayed to its end.
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

Running deviations for a live measurement: one [`StreamingStability`](@ref)
per name in `devs` (`:adev`, `:mdev`, `:hdev`, `:mhdev`), fed through
`push!(live, reading)`. `quantity` says what the readings are:

- `:phase`: time error in seconds;
- `:frequency`: fractional frequency, integrated as `x = [0; cumsum(y)]·τ₀`
  like [`FrequencyData`](@ref);
- `:hz`: counter readings in hertz, turned into fractional frequency
  against `nominal` (by default the first reading to 4 significant
  figures).

`snapshot(live)` returns one `StabilityResult` per deviation; at every
reading count it equals the batch estimators on the readings so far.
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

Number of readings pushed so far.
"""
nsamples(live::LiveStability) = live.readings
