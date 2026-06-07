module SigmaTauTablesExt

using SigmaTau
import Tables

const RESULT_NAMES = (
    :deviation_type,
    :tau,
    :dev,
    :noise_type,
    :ci_lower,
    :ci_upper,
    :edf,
)
const RESULT_TYPES = (
    Symbol,
    Float64,
    Float64,
    Union{Missing,Symbol},
    Union{Missing,Float64},
    Union{Missing,Float64},
    Union{Missing,Float64},
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
    },
}

struct ResultRows
    result::SigmaTau.StabilityResult
end

struct SuiteRows
    suite::SigmaTau.StabilitySuite
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

Base.IteratorSize(::Type{ResultRows}) = Base.HasLength()
Base.length(rows::ResultRows) = length(rows.result.tau)
Base.eltype(::Type{ResultRows}) = RESULT_ROW

Base.IteratorSize(::Type{SuiteRows}) = Base.HasLength()
Base.length(rows::SuiteRows) = sum(result -> length(result.tau), rows.suite.results; init=0)
Base.eltype(::Type{SuiteRows}) = RESULT_ROW

@inline _maybe(values::Vector, i::Int) = i <= length(values) ? values[i] : missing

@inline function _row(result::SigmaTau.StabilityResult, i::Int)::RESULT_ROW
    return (
        deviation_type = result.deviation_type,
        tau = result.tau[i],
        dev = result.dev[i],
        noise_type = _maybe(result.noise_type, i),
        ci_lower = _maybe(result.ci_lower, i),
        ci_upper = _maybe(result.ci_upper, i),
        edf = _maybe(result.edf, i),
    )
end

function Base.iterate(rows::ResultRows, state::Int = 1)
    state > length(rows) && return nothing
    return _row(rows.result, state), state + 1
end

function Base.iterate(rows::SuiteRows, state::Tuple{Int,Int} = (1, 1))
    result_i, row_i = state
    while result_i <= length(rows.suite.results)
        result = rows.suite.results[result_i]
        if row_i <= length(result.tau)
            return _row(result, row_i), (result_i, row_i + 1)
        end
        result_i += 1
        row_i = 1
    end
    return nothing
end

end
