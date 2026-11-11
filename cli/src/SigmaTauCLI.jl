# sigmatau: SigmaTau's command-line interface, one binary with subcommands,
# written to compile with `juliac --trim=safe`. The explorer (app/web and
# app/desktop) runs the same binary with `--format json`.
#
#   julia --project=cli -m SigmaTauCLI analyze clock.dat --devs adev,mdev
#   juliac --output-exe sigmatau --bundle cli/build --trim=safe --experimental cli
#
# Everything stays type-stable for trimming: options are parsed into one
# concrete struct, estimators are chosen through the library's registry (an
# if-chain, not a Dict of functions), and JSON and CSV are written by hand.
# All output is rendered into an IOBuffer first and written in one go.
module SigmaTauCLI

# A program that Julia 1.12 trims cannot start a task, so every threaded loop
# fails in it: refuse to precompile, which a build does first. Pkg skips the
# julia compat when a Manifest is already there. The docs include this file
# on Julia 1.11 without precompiling it.
VERSION < v"1.13" && Base.generating_output() &&
    error("SigmaTauCLI needs Julia 1.13 or later")

using SigmaTau
using Random

# Baked in at compile time; pkgversion needs the source tree at run time.
const SIGMATAU_VERSION = string(pkgversion(SigmaTau))
const SCHEMA_VERSION = 1

# Exit codes.
const EXIT_OK = 0
const EXIT_USAGE = 1
const EXIT_INPUT = 2
const EXIT_ANALYSIS = 3
const EXIT_OUTPUT = 4

struct CLIError <: Exception
    code::Int
    msg::String
end
usage_error(msg::String) = CLIError(EXIT_USAGE, msg)

include("options.jl")
include("record.jl")
include("output.jl")
include("commands.jl")
include("report.jl")
include("live.jl")
include("completion.jl")

# The message for an exception that is not a CLIError. No showerror: it
# dispatches on the exception type, which trimming cannot resolve.
function error_message(err)
    # A threaded loop wraps what its task threw; name that instead.
    if err isa CompositeException && !isempty(err.exceptions)
        err = err.exceptions[1]
        err isa TaskFailedException && (err = err.task.result)
    end
    err isa CLIError && return err.msg
    if err isa ArgumentError
        m = err.msg
        return m isa String ? m : "invalid argument"
    end
    err isa DomainError && return "domain error in the analysis"
    err isa BoundsError && return "index out of range in the analysis"
    return "failed (" * string(nameof(typeof(err))) * ")"
end

# Plain open/close: a do-block closure would box captured locals.
function write_file(path::String, bytes::Vector{UInt8})
    io = open(path, "w")
    try
        write(io, bytes)
    finally
        close(io)
    end
end

function report_error(out::IO, err::IO, o::Options, code::Int, msg::String)
    print(err, "sigmatau: "); print(err, msg); print(err, '\n')
    if o.format == "json"
        buf = IOBuffer()
        write_error_json(buf, o, code, msg)
        write(out, take!(buf))
    end
    return code
end

"""
    sigmatau(args; out, err) -> Int

Run one `sigmatau` command and return its exit code (0 ok, 1 usage, 2 input,
3 analysis, 4 output). Never calls `exit`, so it is safe to call from tests.
"""
function sigmatau(args::Vector{String}; out::IO=stdout, err::IO=stderr)::Int
    if isempty(args) || args[1] in ("help", "-h", "--help")
        if length(args) >= 2 && args[2] != "completion" && !(args[2] in COMMANDS)
            print(err, "sigmatau: no help for \"", args[2], "\"; commands: ", join(COMMANDS, ", "), "\n")
            return EXIT_USAGE
        end
        help_buf = IOBuffer()
        length(args) >= 2 ? print_command_help(help_buf, args[2]) : print_help(help_buf)
        write(out, take!(help_buf))
        return EXIT_OK
    end
    cmd = args[1]
    if cmd in ("version", "--version", "-V")
        print(out, "sigmatau "); print(out, SIGMATAU_VERSION); print(out, '\n')
        return EXIT_OK
    elseif cmd == "completion"
        return run_completion(args, out, err)
    elseif !(cmd in COMMANDS)
        print(err, startswith(cmd, "-") ? "sigmatau: unknown option " : "sigmatau: unknown command ")
        print(err, cmd); print(err, "\nRun 'sigmatau help' for usage.\n")
        return EXIT_USAGE
    end

    o = Options()
    o.command = cmd
    (cmd == "live" || cmd == "model") && (o.format = "json")    # they write JSON only
    cmd == "holdover" && (o.confidence = 0.95)
    for a in args
        a == "--format=json" && (o.format = "json")   # so a usage error can still answer in JSON
    end
    for i in 2:length(args)
        args[i] == "--format" && i < length(args) && args[i + 1] == "json" && (o.format = "json")
    end
    try
        parse_args!(o, args)
    catch e
        return report_error(out, err, o, EXIT_USAGE, error_message(e))
    end
    if o.help
        cmd_buf = IOBuffer()
        print_command_help(cmd_buf, cmd)
        write(out, take!(cmd_buf))
        return EXIT_OK
    end

    if cmd == "live"
        try
            return run_live(out, o)
        catch e
            return report_error(out, err, o, e isa CLIError ? e.code : EXIT_INPUT, error_message(e))
        end
    end

    buf = IOBuffer()
    if cmd == "model"
        try
            run_model(buf, o)
        catch e
            return report_error(out, err, o, e isa CLIError ? e.code : EXIT_ANALYSIS, error_message(e))
        end
    else
        t_start = time_ns()
        data = try
            load_record(o)
        catch e
            return report_error(out, err, o, e isa CLIError ? e.code : EXIT_INPUT, error_message(e))
        end
        try
            run_command(buf, o, data, t_start, args)
        catch e
            return report_error(out, err, o, e isa CLIError ? e.code : EXIT_ANALYSIS, error_message(e))
        end
    end
    bytes = take!(buf)
    try
        if isempty(o.out)
            write(out, bytes)
        else
            write_file(o.out, bytes)
            o.quiet || (print(err, "sigmatau: wrote "); print(err, o.out); print(err, '\n'))
        end
    catch e
        # A consumer that closed the pipe early (`sigmatau analyze … | head`) is not a failure.
        e isa Base.IOError && isempty(o.out) && return EXIT_OK
        return report_error(out, err, o, EXIT_OUTPUT, "could not write the output: " * error_message(e))
    end
    return EXIT_OK
end

# Standard output as file descriptor 1. A write to Core.stdout goes to libuv,
# which reports a closed pipe on stderr (`jl_uv_writecb() ERROR: broken pipe`)
# and throws nothing, so `live` never saw its reader go. A write here fails
# with an IOError, which the commands already treat as the reader's goodbye.
struct FdOut <: IO
    fd::Cint
end

function Base.unsafe_write(io::FdOut, p::Ptr{UInt8}, n::UInt)
    done = 0
    while done < n
        w = ccall(:write, Cssize_t, (Cint, Ptr{UInt8}, Csize_t), io.fd, p + done, n - done)
        if w >= 0
            done += w
            continue
        end
        errno = Libc.errno()
        errno == Libc.EINTR && continue
        # Another part of the runtime may have made the pipe non-blocking: wait for room.
        errno == Libc.EAGAIN && (Libc.systemsleep(0.001); continue)
        throw(Base.IOError("write", -errno))
    end
    return Int(n)
end

function Base.write(io::FdOut, b::UInt8)
    r = Ref(b)
    GC.@preserve r unsafe_write(io, Base.unsafe_convert(Ptr{UInt8}, r), UInt(1))
    return 1
end

Base.flush(::FdOut) = nothing

function (@main)(args::Vector{String})::Cint
    # Windows has no write(2) in the runtime: keep Core.stdout there.
    Sys.iswindows() && return Cint(sigmatau(args; out=Core.stdout, err=Core.stderr))
    return Cint(sigmatau(args; out=FdOut(Cint(1)), err=Core.stderr))
end

end # module
