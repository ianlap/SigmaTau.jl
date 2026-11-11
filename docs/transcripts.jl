# Command-line transcripts for the docs pages. `sigmatau("analyze clock.dat")`
# prints the command line and what the program answers, with the numbers cut
# to four significant digits. `sigmatau_json` prints a JSON answer indented,
# with long arrays cut short, whole or one member of it. The commands run in
# process, in a scratch directory that holds the sample files.

using SigmaTau, Printf, Random

include(joinpath(@__DIR__, "..", "cli", "src", "SigmaTauCLI.jl"))

const CLI_DIR = mktempdir()

function write_columns(name, columns...)
    open(joinpath(CLI_DIR, name), "w") do io
        for row in zip(columns...)
            join(io, row, ' ')
            println(io)
        end
    end
end

# clock.dat: time and phase at 1 s, with a frequency offset and ten samples
# missing. a.dat and b.dat: one device against two references, in frequency.
let N = 4096, t = 0:N-1, keep = setdiff(1:N, 1000:1009)
    p = noise_gen(PhaseData, N, 1.0; sigma1 = Dict(0 => 1e-12, -2 => 3e-14),
                  rng = Xoshiro(5))
    write_columns("clock.dat", t[keep], (p.x .+ 1e-11 .* t)[keep])
    wfm(h, seed) = noise_gen(FrequencyData, N, 1.0; h = Dict(:WHFM => h),
                             rng = Xoshiro(seed)).y
    device = wfm(1e-24, 2)
    write_columns("a.dat", t, device .- wfm(1e-23, 3))
    write_columns("b.dat", t, device .- wfm(1e-23, 4))
end

# ocxo.txt: the phase noise of a 10 MHz oscillator, offset in Hz and ℒ(f)
# in dBc/Hz.
write_columns("ocxo.txt", [1, 10, 100, 1e3, 1e4, 1e5], [-100, -130, -150, -160, -165, -165])

# A number with a decimal point that stands as a CSV field or a JSON value.
const NUMBER = r"(?<=[,:\[\n])-?\d+\.\d+(?:e-?\d+)?(?=[,\]}\n])"

function answer(line)
    io = IOBuffer()
    cd(() -> SigmaTauCLI.sigmatau(String.(split(line)); out = io, err = io), CLI_DIR)
    short(x) = string(parse(Float64, @sprintf("%.3e", parse(Float64, x))))
    return replace(String(take!(io)), NUMBER => short)
end

# Break a long line at a blank, as a terminal does.
function wrapped(line; width = 88)
    length(line) <= width && return line
    i = findprev(' ', line, thisind(line, width))
    i === nothing && return line
    return line[1:i-1] * '\n' * wrapped(line[i+1:end]; width)
end

sigmatau(line) = print("\$ sigmatau ", line, '\n',
                       join(wrapped.(split(answer(line), '\n')), '\n'))

# Indent a JSON text. An array of scalars stays on one line, cut to its first
# `keep` elements, and an array of arrays or objects to its first element.
function indent_json(json; keep = 3)
    c, out, depth, i = collect(json), IOBuffer(), 0, 1
    newline() = print(out, '\n', "  "^depth)
    while i <= length(c)
        if c[i] == '"'
            j = findnext(k -> c[k] == '"' && c[k-1] != '\\', eachindex(c), i + 1)
            print(out, String(c[i:j])); i = j
        elseif c[i] == '[' && c[i+1] in "[{"
            j, k = closing(c, i), closing(c, i + 1)
            depth += 1; print(out, '['); newline()
            print(out, replace(indent_json(String(c[i+1:k]); keep), "\n" => "\n" * "  "^depth))
            k + 1 < j && (print(out, ','); newline(); print(out, '…'))
            depth -= 1; newline(); print(out, ']'); i = j
        elseif c[i] == '['
            j = closing(c, i)
            items = split(String(c[i+1:j-1]), ','; keepempty = false)
            print(out, '[', join(first(items, keep), ", "), length(items) > keep ? ", …]" : "]"); i = j
        elseif c[i] == '{'
            depth += 1; print(out, '{'); newline()
        elseif c[i] == '}'
            depth -= 1; newline(); print(out, '}')
        elseif c[i] == ','
            print(out, ','); newline()
        elseif c[i] == ':'
            print(out, ": ")
        else
            print(out, c[i])
        end
        i += 1
    end
    return String(take!(out))
end

# Index of the bracket that closes the one at `i`.
function closing(c, i)
    depth = 0
    return findnext(k -> (depth += (c[k] in "[{") - (c[k] in "]}")) == 0, eachindex(c), i)
end

# With `only`, that member of the answer alone, the others shown as `…`.
function sigmatau_json(line; only = "", keep = 3)
    text = strip(answer(line))
    json = text[findfirst('{', text):end]      # an error line comes first
    print("\$ sigmatau ", line, '\n', chopsuffix(text, json))
    isempty(only) && return println(indent_json(json; keep))
    c = collect(json)
    key = collect("\"$only\":")
    i = findfirst(k -> c[k:k+length(key)-1] == key, 1:length(c)-length(key)) + length(key)
    value = indent_json(String(c[i:closing(c, i)]); keep)
    println("{\n  …,\n  \"", only, "\": ", replace(value, "\n" => "\n  "), ",\n  …\n}")
end
