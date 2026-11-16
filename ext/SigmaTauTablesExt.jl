module SigmaTauTablesExt

using SigmaTau
import Tables

# `ci_lower` / `ci_upper` stay flat numeric columns (the per-row split of the
# `ci` tuple field, matching the `ci_lower` / `ci_upper` accessors); `neff` is
# always populated.
const RESULT_NAMES = (
    :deviation_type,
    :tau,
    :dev,
    :noise_type,
    :ci_lower,
    :ci_upper,
    :edf,
    :neff,
)
const RESULT_TYPES = (
    Symbol,
    Float64,
    Float64,
    Union{Missing,Symbol},
    Union{Missing,Float64},
    Union{Missing,Float64},
    Union{Missing,Float64},
    Int,
)
const RESULT_SCHEMA = Tables.Schema(RESULT_NAMES, RESULT_TYPES)
const RESULT_ROW = NamedTuple{
    RESULT_NAMES,
    Tuple{
        Symbol,
        Float64,
        Float64,
        Union{Missing,Symbol},
        Union{Missing,Float64},
        Union{Missing,Float64},
        Union{Missing,Float64},
        Int,
    },
}

struct ResultRows
    result::SigmaTau.StabilityResult
end

struct SuiteRows
    suite::SigmaTau.StabilitySuite
end

# A hat's rows: the clock's name, then the columns of a result.
const HAT_NAMES = (:clock, RESULT_NAMES...)
const HAT_SCHEMA = Tables.Schema(HAT_NAMES, (String, RESULT_TYPES...))
const HAT_ROW = NamedTuple{HAT_NAMES, Tuple{String, fieldtypes(RESULT_ROW)...}}

struct HatRows
    hat::SigmaTau.HatResult
end

Tables.istable(::Type{SigmaTau.StabilityResult}) = true
Tables.rowaccess(::Type{SigmaTau.StabilityResult}) = true
Tables.rows(result::SigmaTau.StabilityResult) = ResultRows(result)
Tables.schema(::SigmaTau.StabilityResult) = RESULT_SCHEMA
Tables.schema(::ResultRows) = RESULT_SCHEMA

Tables.istable(::Type{SigmaTau.StabilitySuite}) = true
Tables.rowaccess(::Type{SigmaTau.StabilitySuite}) = true
Tables.rows(suite::SigmaTau.StabilitySuite) = SuiteRows(suite)
Tables.schema(::SigmaTau.StabilitySuite) = RESULT_SCHEMA
Tables.schema(::SuiteRows) = RESULT_SCHEMA

Tables.istable(::Type{SigmaTau.HatResult}) = true
Tables.rowaccess(::Type{SigmaTau.HatResult}) = true
Tables.rows(hat::SigmaTau.HatResult) = HatRows(hat)
Tables.schema(::SigmaTau.HatResult) = HAT_SCHEMA
Tables.schema(::HatRows) = HAT_SCHEMA

Base.IteratorSize(::Type{ResultRows}) = Base.HasLength()
Base.length(rows::ResultRows) = length(rows.result.tau)
Base.eltype(::Type{ResultRows}) = RESULT_ROW

Base.IteratorSize(::Type{SuiteRows}) = Base.HasLength()
Base.length(rows::SuiteRows) = sum(result -> length(result.tau), rows.suite.results; init=0)
Base.eltype(::Type{SuiteRows}) = RESULT_ROW

Base.IteratorSize(::Type{HatRows}) = Base.HasLength()
Base.length(rows::HatRows) = sum(result -> length(result.tau), rows.hat.clocks; init=0)
Base.eltype(::Type{HatRows}) = HAT_ROW

@inline _maybe(values::Vector, i::Int) = i <= length(values) ? values[i] : missing

@inline function _row(result::SigmaTau.StabilityResult, i::Int)::RESULT_ROW
    ci = i <= length(result.ci) ? result.ci[i] : missing
    return (
        deviation_type = result.deviation_type,
        tau = result.tau[i],
        dev = result.dev[i],
        noise_type = _maybe(result.noise_type, i),
        ci_lower = ci === missing ? missing : ci.lo,
        ci_upper = ci === missing ? missing : ci.hi,
        edf = _maybe(result.edf, i),
        neff = result.neff[i],
    )
end

function Base.iterate(rows::ResultRows, state::Int = 1)
    state > length(rows) && return nothing
    return _row(rows.result, state), state + 1
end

# Row `row_i` of result `result_i` in a sequence of results, or the next
# row there is; `nothing` past the last.
function _next(results, (result_i, row_i)::Tuple{Int,Int})
    while result_i <= length(results)
        row_i <= length(results[result_i].tau) && return result_i, row_i
        result_i += 1
        row_i = 1
    end
    return nothing
end

function Base.iterate(rows::SuiteRows, state::Tuple{Int,Int} = (1, 1))
    at = _next(rows.suite.results, state)
    at === nothing && return nothing
    result_i, row_i = at
    return _row(rows.suite.results[result_i], row_i), (result_i, row_i + 1)
end

function Base.iterate(rows::HatRows, state::Tuple{Int,Int} = (1, 1))
    at = _next(rows.hat.clocks, state)
    at === nothing && return nothing
    clock_i, row_i = at
    row = HAT_ROW((rows.hat.names[clock_i], _row(rows.hat.clocks[clock_i], row_i)...))
    return row, (clock_i, row_i + 1)
end

end
