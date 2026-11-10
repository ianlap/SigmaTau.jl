# Byte-level reading of large text records: a chunked line reader, field
# splitting and number parsing on byte ranges, so a file of any length costs
# its output arrays and one read buffer, not a String per line or field.

mutable struct _Lines
    io::IOStream
    buf::Vector{UInt8}      # buf[lo:hi] is read but not yet returned
    chunk::Vector{UInt8}
    lo::Int
    hi::Int
    done::Bool              # the stream is exhausted
end

_Lines(path::AbstractString) =
    _Lines(open(path, "r"), Vector{UInt8}(undef, 1 << 20), Vector{UInt8}(undef, 1 << 20), 1, 0, false)

Base.close(l::_Lines) = close(l.io)

function _refill!(l::_Lines)
    rest = l.hi - l.lo + 1
    rest > 0 && copyto!(l.buf, 1, l.buf, l.lo, rest)
    n = readbytes!(l.io, l.chunk)
    n == 0 && (l.done = true)
    length(l.buf) < rest + n && resize!(l.buf, 2 * (rest + n))
    copyto!(l.buf, rest + 1, l.chunk, 1, n)
    l.lo = 1; l.hi = rest + n
    return l
end

# The next line as a byte range of l.buf, without its "\n" or "\r\n";
# `nothing` at the end of the file.
function _nextline!(l::_Lines)::Union{UnitRange{Int},Nothing}
    while true
        nl = 0
        @inbounds for k in l.lo:l.hi
            if l.buf[k] == 0x0a
                nl = k
                break
            end
        end
        if nl > 0 || (l.done && l.lo <= l.hi)
            s = l.lo
            e = nl > 0 ? nl - 1 : l.hi
            l.lo = e + 2
            (e >= s && l.buf[e] == 0x0d) && (e -= 1)
            return s:e
        end
        l.done && return nothing
        _refill!(l)
    end
end

_isblank(b::UInt8) = b == 0x20 || b == 0x09 || b == 0x0d

function _strip(buf::Vector{UInt8}, r::UnitRange{Int})
    s, e = first(r), last(r)
    while s <= e && _isblank(buf[s]); s += 1; end
    while e >= s && _isblank(buf[e]); e -= 1; end
    # A double-quoted field ("1.0").
    if e - s >= 1 && buf[s] == 0x22 && buf[e] == 0x22
        s += 1; e -= 1
    end
    return s:e
end

_blank_or_comment(buf::Vector{UInt8}, r::UnitRange{Int}) =
    (s = _strip(buf, r); isempty(s) || buf[first(s)] == UInt8('#'))

# Split the line r of buf into fields, written to `out` (emptied first): at
# each `delim`, or at runs of blanks when `delim` is '\0'. Fields are stripped.
function _fields!(out::Vector{UnitRange{Int}}, buf::Vector{UInt8}, r::UnitRange{Int}, delim::Char)
    empty!(out)
    if delim == '\0'
        k = first(r)
        while k <= last(r)
            while k <= last(r) && _isblank(buf[k]); k += 1; end
            k > last(r) && break
            s = k
            while k <= last(r) && !_isblank(buf[k]); k += 1; end
            push!(out, _strip(buf, s:k-1))
        end
    else
        d = UInt8(delim)
        s = first(r)
        for k in r
            if buf[k] == d
                push!(out, _strip(buf, s:k-1))
                s = k + 1
            end
        end
        push!(out, _strip(buf, s:last(r)))
    end
    return out
end

# The number in buf[r], or `nothing`. The C library's strtod, as
# `parse(Float64, ::String)` uses, on the bytes in place.
function _parse_f64(buf::Vector{UInt8}, r::UnitRange{Int})::Union{Float64,Nothing}
    isempty(r) && return nothing
    ok, v = GC.@preserve buf ccall(:jl_try_substrtod, Tuple{Bool,Float64},
                                   (Ptr{UInt8}, Csize_t, Csize_t), pointer(buf), first(r) - 1, length(r))
    return ok ? v : nothing
end

_field_string(buf::Vector{UInt8}, r::UnitRange{Int}) = String(buf[r])

function _digits(buf::Vector{UInt8}, s::Int, n::Int)
    v = 0
    for k in s:s+n-1
        b = buf[k]
        0x30 <= b <= 0x39 || return -1
        v = 10v + (b - 0x30)
    end
    return v
end

# UTC seconds since 2000-01-01 of an ISO 8601 timestamp in buf[r]
# ("2026-11-08T21:41:03.5", a space for the T, then Z, ±hh:mm, ±hhmm or
# nothing); NaN if it is not one.
function _iso_seconds(buf::Vector{UInt8}, r::UnitRange{Int})
    s = first(r); n = length(r)
    n >= 19 || return NaN
    (buf[s+4] == UInt8('-') && buf[s+7] == UInt8('-') && buf[s+13] == UInt8(':') &&
     buf[s+16] == UInt8(':') && (buf[s+10] == UInt8('T') || buf[s+10] == UInt8(' '))) || return NaN
    y = _digits(buf, s, 4); mo = _digits(buf, s + 5, 2); d = _digits(buf, s + 8, 2)
    h = _digits(buf, s + 11, 2); mi = _digits(buf, s + 14, 2); se = _digits(buf, s + 17, 2)
    min(y, mo, d, h, mi, se) < 0 && return NaN
    (1 <= mo <= 12 && 1 <= d <= 31) || return NaN
    k = s + 19
    frac = 0.0
    if k <= last(r) && buf[k] == UInt8('.')
        e = k + 1
        while e <= last(r) && 0x30 <= buf[e] <= 0x39; e += 1; end
        f = _parse_f64(buf, k:e-1)
        f === nothing && return NaN
        frac = f; k = e
    end
    offset = 0
    if k <= last(r)
        b = buf[k]
        if b == UInt8('Z')
            k += 1
        elseif b == UInt8('+') || b == UInt8('-')
            hh = k + 2 <= last(r) ? _digits(buf, k + 1, 2) : -1
            j = k + 3
            j <= last(r) && buf[j] == UInt8(':') && (j += 1)
            mm = j + 1 <= last(r) ? _digits(buf, j, 2) : -1
            min(hh, mm) < 0 && return NaN
            offset = (b == UInt8('+') ? 1 : -1) * (hh * 3600 + mm * 60)
            k = j + 2
        end
    end
    k == last(r) + 1 || return NaN
    days = Dates.value(Dates.Date(y, mo, d) - Dates.Date(2000, 1, 1))
    return days * 86400.0 + h * 3600 + mi * 60 + se + frac - offset
end

_iso_seconds(s::AbstractString) = (b = Vector{UInt8}(codeunits(strip(s))); _iso_seconds(b, 1:length(b)))
