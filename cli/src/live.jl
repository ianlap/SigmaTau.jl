# live: readings from an instrument or the simulator turned into running
# σ(τ), written as one JSON object per line while the measurement runs.

# The τ grid when --count is not given: m up to 2^16, a ring buffer of
# 4·2^16 samples per deviation.
const LIVE_GRID_N = 1 << 17

# The simulated clock: white PM at short τ, white FM in the middle, random
# walk FM at long τ (levels per sample; see SimulatedSource).
simulator(o::Options, tau0::Float64) =
    SimulatedSource(white_pm=1e-11, white_fm=2e-12, rw_fm=3e-16, tau0=tau0, rng=Xoshiro(o.seed))

function run_live(out::IO, o::Options)::Int
    for d in o.devs
        haskey(SigmaTau._STREAMING_LAG, d) || throw(usage_error("live streams adev, mdev, hdev and mhdev, not " * String(d)))
    end
    tau0 = isnan(o.tau0) ? 1.0 : o.tau0
    ms = isempty(o.m) ? tau_values(tau_mode(o.grid), o.count > 0 ? max(o.count, 4) : LIVE_GRID_N, :adev) : o.m
    quantity = o.hz ? :hz : o.freq ? :frequency : :phase
    live = LiveStability(o.devs, tau0, ms; quantity=quantity, nominal=o.nominal)
    src = o.input
    if src == "sim"
        quantity === :phase || throw(usage_error("the simulator reads phase; drop --freq and --hz"))
        return live_loop(out, o, live, simulator(o, tau0), !o.fast)
    elseif startswith(src, "tcp://")
        Sys.iswindows() && throw(usage_error("tcp:// sources need macOS or Linux; on Windows use a serial port"))
        host, port = host_port(src)
        sock = tcp_connect(host, port)
        try
            return live_loop(out, o, live, InstrumentSource(sock, o.query == "\0" ? "READ?" : o.query), false)
        finally
            close(sock)
        end
    end
    query = o.query == "\0" ? "" : o.query
    ispath(src) || throw(CLIError(EXIT_INPUT, "no such source: " * src))
    io = isempty(query) ? open(src, "r") : Duplex(open(src, "r"), open(src; write=true, create=false, truncate=false))
    try
        return live_loop(out, o, live, InstrumentSource(io, query), false)
    finally
        close(io)
    end
end

function host_port(src::String)
    parts = rsplit(src[7:end], ':'; limit=2)
    port = length(parts) == 2 ? tryparse(Int, parts[2]) : nothing
    (port === nothing || !(0 < port < 65536) || isempty(parts[1])) &&
        throw(usage_error("give the source as tcp://HOST:PORT, got " * src))
    return String(parts[1]), port
end

# A two-way stream as two IOStreams: one IOStream that is both read and
# written sends its read buffer back out with the next write.
struct Duplex <: IO
    r::IOStream
    w::IOStream
end
Base.eof(d::Duplex) = eof(d.r)
Base.readline(d::Duplex; keep::Bool=false) = readline(d.r; keep=keep)
Base.write(d::Duplex, b::UInt8) = write(d.w, b)
Base.unsafe_write(d::Duplex, p::Ptr{UInt8}, n::UInt) = unsafe_write(d.w, p, n)
Base.flush(d::Duplex) = flush(d.w)
Base.close(d::Duplex) = (close(d.r); close(d.w))

# An IPv4 TCP connection through libc: the Sockets library reads through
# libuv callbacks, which trimming cannot compile.
function tcp_connect(host::String, port::Int)::Duplex
    fail(why::String) = CLIError(EXIT_INPUT, "could not connect to " * host * ":" * string(port) * ": " * why)
    he = ccall(:gethostbyname, Ptr{UInt8}, (Cstring,), host)
    he == C_NULL && throw(fail("unknown host"))
    addr = unsafe_load(unsafe_load(Ptr{Ptr{Ptr{UInt8}}}(he + 24)))   # hostent.h_addr_list[0]
    addr == C_NULL && throw(fail("no IPv4 address"))
    sa = zeros(UInt8, 16)                       # sockaddr_in, AF_INET = 2
    Sys.isapple() || Sys.isbsd() ? (sa[1] = 16; sa[2] = 2) : (sa[1] = 2)
    sa[3] = UInt8(port >> 8); sa[4] = UInt8(port & 0xff)
    for k in 1:4
        sa[4 + k] = unsafe_load(addr, k)
    end
    fd = ccall(:socket, Cint, (Cint, Cint, Cint), 2, 1, 0)   # SOCK_STREAM = 1
    fd < 0 && throw(fail("no socket"))
    if ccall(:connect, Cint, (Cint, Ptr{UInt8}, UInt32), fd, sa, 16) != 0
        ccall(:close, Cint, (Cint,), fd)
        throw(fail("connection refused or unreachable"))
    end
    return Duplex(fdio(fd, true), fdio(ccall(:dup, Cint, (Cint,), fd), true))
end

# Read until --count or the end of the stream, writing a line every --every
# readings and a last one with "done": true. A reader that goes away ends
# the run quietly.
function live_loop(out::IO, o::Options, live::LiveStability, src, paced::Bool)::Int
    t_start = time_ns()
    batch = Float64[]
    while o.count == 0 || nsamples(live) < o.count
        if paced
            wait_s = nsamples(live) * live.tau0 - (time_ns() - t_start) / 1e9
            wait_s > 0 && Libc.systemsleep(wait_s)
        end
        r = read_reading(src)
        isnan(r) && break
        push!(live, r)
        push!(batch, r)
        if length(batch) >= o.every && nsamples(live) != o.count
            live_line(out, o, live, batch, t_start, false) || return EXIT_OK
            empty!(batch)
        end
    end
    live_line(out, o, live, batch, t_start, true)
    return EXIT_OK
end

function live_line(out::IO, o::Options, live::LiveStability, batch::Vector{Float64},
                   t_start::UInt64, done::Bool)::Bool
    io = IOBuffer()
    json_head(io, o, o.input)
    json_key(io, "quantity"); json_string(io, String(live.quantity))
    json_key(io, "tau0"); json_number(io, live.tau0)
    json_key(io, "nominal"); json_number(io, live.nominal)
    json_key(io, "n"); json_number(io, nsamples(live))
    json_key(io, "readings"); json_array(io, batch)
    results = snapshot(live)
    json_key(io, "tau"); json_array(io, results[1].tau)
    json_key(io, "results"); print(io, '[')
    for (k, r) in enumerate(results)
        k > 1 && print(io, ',')
        print(io, '{')
        json_key(io, "deviation", true); json_string(io, String(r.deviation_type))
        json_key(io, "dev"); json_array(io, r.dev)
        json_key(io, "neff"); json_array(io, r.neff)
        print(io, '}')
    end
    print(io, ']')
    json_key(io, "done"); json_bool(io, done)
    json_key(io, "elapsed_ms"); json_number(io, (time_ns() - t_start) / 1e6)
    print(io, "}\n")
    try
        write(out, take!(io))
        flush(out)
    catch e
        e isa Base.IOError && return false
        rethrow()
    end
    return true
end
