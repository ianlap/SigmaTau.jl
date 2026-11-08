# Names replaced in 0.6. Each forwards to its replacement with a deprecation
# warning (shown with `julia --depwarn=yes`, which `Pkg.test` turns on) and
# will be removed in 0.7.

function Sy(data::AbstractTimingData; window::Symbol = :hann, taper::Symbol = window, kwargs...)
    Base.depwarn("`Sy` is deprecated, use `psd_y` (with `taper=` for `window=`)", :Sy)
    return psd_y(data; taper = taper, kwargs...)
end

function Sx(data::AbstractTimingData; window::Symbol = :hann, taper::Symbol = window, kwargs...)
    Base.depwarn("`Sx` is deprecated, use `psd_x` (with `taper=` for `window=`)", :Sx)
    return psd_x(data; taper = taper, kwargs...)
end

function L(data::AbstractTimingData; window::Symbol = :hann, taper::Symbol = window, kwargs...)
    Base.depwarn("`L` is deprecated, use `ssb_phase_noise` (with `taper=` for `window=`)", :L)
    return ssb_phase_noise(data; taper = taper, kwargs...)
end
export Sy, Sx, L

Base.@deprecate nch(args...) ncornered_hat(args...)

Base.@deprecate save_result(path::AbstractString, r::StabilityResult) save(path, r)
Base.@deprecate save_suite(path::AbstractString, s::StabilitySuite; kwargs...) save(path, s; kwargs...)
Base.@deprecate load_result(path::AbstractString) load(path)
Base.@deprecate load_suite(path::AbstractString) load(path)

function Base.getproperty(r::SpectralResult, s::Symbol)
    if s === :window
        Base.depwarn("`SpectralResult.window` is deprecated, use `.taper`", :getproperty)
        return getfield(r, :taper)
    end
    return getfield(r, s)
end
