module MQLib

using Printf

import MQLib_jll
import QUBOTools
import QUBODrivers
import MathOptInterface as MOI

const __VERSION__ = v"0.1.0"
const _HEURISTICS = Dict{String,String}()
const _MQLIB_C_ABI_VERSION = UInt32(1)
const _MQLIB_C_INDEX_BASE_ONE = Int32(1)
const _MQLIB_STATUS_OK = Cint(0)
const _MQLIB_STATUS_BUFFER_TOO_SMALL = Cint(4)

struct _MQLibCQUBOInput
    abi_version::UInt32
    dimension::Int32
    linear::Ptr{Cdouble}
    quadratic_count::Int64
    quadratic_i::Ptr{Int32}
    quadratic_j::Ptr{Int32}
    quadratic_value::Ptr{Cdouble}
    index_base::Int32
    heuristic::Cstring
    runtime_limit_seconds::Cdouble
    random_seed::Int32
    hyperheuristic_data_dir::Cstring
end

struct _MQLibCQUBOResult
    abi_version::UInt32
    objective_value::Cdouble
    runtime_seconds::Cdouble
    solution::Ptr{Int32}
    solution_length::Int32
    selected_heuristic::Cstring
    selected_heuristic_length::Int32
    history_objective_values::Ptr{Cdouble}
    history_times_seconds::Ptr{Cdouble}
    history_capacity::Int32
    history_length::Int32
end

function __init__()
    let exe = MQLib_jll.MQLib()
        ms = eachmatch(r"([a-zA-Z0-9]+)\r?\n\s+([^\r\n]+)\r?\n?", read(`$exe -l`, String))

        for m in ms
            push!(_HEURISTICS, m[1] => m[2])
        end
    end

    return nothing
end

QUBODrivers.@setup Optimizer begin
    name       = "MQLib"
    version    = __VERSION__
    attributes = begin
        RandomSeed["seed"]::Union{Integer,Nothing} = nothing
        NumberOfReads["num_reads"]::Integer = 1
        Heuristic["heuristic"]::Union{String,Nothing} = nothing
    end
end

function QUBODrivers.sample(sampler::Optimizer{T}) where {T}
    n, L, Q, α, β = QUBOTools.qubo(sampler, :dict; sense = :max, domain = :bool)

    num_reads      = MOI.get(sampler, MQLib.NumberOfReads())
    silent         = MOI.get(sampler, MOI.Silent())
    heuristic      = MOI.get(sampler, MQLib.Heuristic())
    random_seed    = MOI.get(sampler, MQLib.RandomSeed())
    time_limit_sec = MOI.get(sampler, MOI.TimeLimitSec())

    if num_reads <= 0
        error("Number of reads must be a positive integer")
    end

    if !isnothing(heuristic) && !haskey(_HEURISTICS, heuristic)
        error("Invalid QUBO Heuristic code '$heuristic'")
    end

    if !isnothing(random_seed)
        random_seed %= 65_536
    end

    run_time_limit = if isnothing(time_limit_sec)
        1.0 / num_reads
    else
        time_limit_sec / num_reads
    end

    metadata = Dict{String,Any}(
        "time"   => Dict{String,Any}(),
        "origin" => Dict{String,Any}(
            "name"      => "MQLib",
            "version"   => __VERSION__,
            "heuristic" => heuristic,
        ),
    )

    samples, effective_time = if _mqlib_has_c_api()
        _sample_with_c_api(
            T,
            n,
            L,
            Q,
            α,
            β;
            num_reads,
            silent,
            heuristic,
            random_seed,
            run_time_limit,
        )
    else
        _sample_with_executable(
            T,
            n,
            L,
            Q,
            α,
            β;
            num_reads,
            silent,
            heuristic,
            random_seed,
            run_time_limit,
        )
    end

    metadata["time"]["effective"] = effective_time

    return QUBOTools.SampleSet{T}(samples, metadata; sense = :max, domain = :bool)
end

function _sample_with_executable(
    ::Type{T},
    n::Integer,
    L,
    Q,
    α,
    β;
    num_reads::Integer,
    silent::Bool,
    heuristic::Union{String,Nothing},
    random_seed::Union{Integer,Nothing},
    run_time_limit::Float64,
) where {T}
    V = Set{Int}(1:n)
    model = QUBOTools.Model{Int,T,Int}(
        V, L, Q;
        scale = α,
        offset = β,
        sense = :max,
        domain = :bool,
    )

    samples = QUBODrivers.Sample{T,Int}[]
    effective_time = mktempdir() do temp_path
        file_path = joinpath(temp_path, "model.qubo")

        args = _mqlib_args(;
            file_path,
            heuristic,
            random_seed,
            run_time_limit,
        )

        QUBOTools.write_model(file_path, model, QUBOTools.Format{:qubo}(; style = :mqlib))

        let exe = MQLib_jll.MQLib()
            cmd = `$exe $args`

            _print_header(silent, heuristic)

            t = 0.0

            for i = 1:num_reads
                lines = readlines(cmd)
                info  = split(lines[begin], ',')

                λ = parse(T, info[4])
                ψ = parse.(Int, split(lines[end], ' '))
                s = QUBODrivers.Sample{T}(ψ, α * (λ + β))

                push!(samples, s)

                m = collect(eachmatch(r"(([0-9]+):([0-9]+))+", info[6]))
                λ̄ = parse.(Float64, getindex.(m, 2))
                t̄ = parse.(Float64, getindex.(m, 3))
                t += parse(Float64, info[5])

                _print_iter(silent, i, λ̄, t .+ t̄)
            end

            _print_footer(silent)

            t
        end
    end

    return samples, effective_time
end

function _sample_with_c_api(
    ::Type{T},
    n::Integer,
    L,
    Q,
    α,
    β;
    num_reads::Integer,
    silent::Bool,
    heuristic::Union{String,Nothing},
    random_seed::Union{Integer,Nothing},
    run_time_limit::Float64,
) where {T}
    linear, quadratic_i, quadratic_j, quadratic_value = _mqlib_problem_data(n, L, Q)
    samples = QUBODrivers.Sample{T,Int}[]

    _print_header(silent, heuristic)

    t = 0.0
    for i = 1:num_reads
        result = _mqlib_solve_qubo(
            n,
            linear,
            quadratic_i,
            quadratic_j,
            quadratic_value;
            heuristic,
            random_seed = _mqlib_run_seed(random_seed),
            run_time_limit,
        )

        λ = T(result.objective_value)
        ψ = Int.(result.solution)
        push!(samples, QUBODrivers.Sample{T}(ψ, α * (λ + β)))

        _print_iter(
            silent,
            i,
            result.history_objective_values,
            t .+ result.history_times_seconds,
        )
        t += result.runtime_seconds
    end

    _print_footer(silent)

    return samples, t
end

function _mqlib_has_c_api()
    return isdefined(MQLib_jll, :libmqlib_c_api)
end

function _mqlib_library()
    return getproperty(MQLib_jll, :libmqlib_c_api)
end

function _mqlib_hyperheuristic_data_dir()
    return joinpath(MQLib_jll.artifact_dir, "share", "mqlib", "hhdata")
end

function _mqlib_problem_data(n::Integer, L, Q)
    linear = zeros(Cdouble, n)

    for (i, value) in L
        index = Int(i)
        1 <= index <= n || error("Invalid QUBO linear index '$i'")
        linear[index] += Cdouble(value)
    end

    quadratic_i = Int32[]
    quadratic_j = Int32[]
    quadratic_value = Cdouble[]

    for (indices, value) in Q
        i, j = Tuple(indices)
        first = Int(i)
        second = Int(j)
        if first == second
            1 <= first <= n || error("Invalid QUBO quadratic index '$i'")
            linear[first] += Cdouble(value)
        else
            1 <= first <= n || error("Invalid QUBO quadratic index '$i'")
            1 <= second <= n || error("Invalid QUBO quadratic index '$j'")
            if second < first
                first, second = second, first
            end
            push!(quadratic_i, Int32(first))
            push!(quadratic_j, Int32(second))
            push!(quadratic_value, Cdouble(value))
        end
    end

    return linear, quadratic_i, quadratic_j, quadratic_value
end

function _mqlib_run_seed(random_seed::Union{Integer,Nothing})
    if isnothing(random_seed)
        return Int32(rand(0:65_535))
    end

    return Int32(mod(random_seed, 65_536))
end

function _mqlib_cstring(value::Union{AbstractString,Nothing})
    bytes = Vector{UInt8}(codeunits(something(value, "")))
    push!(bytes, 0x00)
    return bytes
end

function _mqlib_status_message(status::Integer)
    message = ccall(
        (:mqlib_c_status_message, _mqlib_library()),
        Cstring,
        (Cint,),
        Cint(status),
    )

    return unsafe_string(message)
end

function _mqlib_solve_qubo(
    n::Integer,
    linear::Vector{Cdouble},
    quadratic_i::Vector{Int32},
    quadratic_j::Vector{Int32},
    quadratic_value::Vector{Cdouble};
    heuristic::Union{String,Nothing},
    random_seed::Integer,
    run_time_limit::Float64,
)
    selected_capacity = Int32(64)
    history_capacity = Int32(1024)

    while true
        call = _mqlib_call_solve_qubo(
            n,
            linear,
            quadratic_i,
            quadratic_j,
            quadratic_value;
            heuristic,
            random_seed,
            run_time_limit,
            selected_capacity,
            history_capacity,
        )

        status = call.status
        result = call.result
        if status == _MQLIB_STATUS_BUFFER_TOO_SMALL &&
           (result.selected_heuristic_length > selected_capacity ||
            result.history_length > history_capacity)
            selected_capacity = max(selected_capacity, result.selected_heuristic_length)
            history_capacity = max(history_capacity, result.history_length)
            continue
        elseif status != _MQLIB_STATUS_OK
            error("MQLib C API failed: $(_mqlib_status_message(status))")
        end

        history_length = Int(result.history_length)
        return (
            objective_value = result.objective_value,
            runtime_seconds = result.runtime_seconds,
            solution = copy(call.solution[1:Int(result.solution_length)]),
            selected_heuristic = _mqlib_selected_heuristic(
                call.selected_heuristic,
                result.selected_heuristic_length,
            ),
            history_objective_values = copy(call.history_objective_values[1:history_length]),
            history_times_seconds = copy(call.history_times_seconds[1:history_length]),
        )
    end
end

function _mqlib_call_solve_qubo(
    n::Integer,
    linear::Vector{Cdouble},
    quadratic_i::Vector{Int32},
    quadratic_j::Vector{Int32},
    quadratic_value::Vector{Cdouble};
    heuristic::Union{String,Nothing},
    random_seed::Integer,
    run_time_limit::Float64,
    selected_capacity::Integer,
    history_capacity::Integer,
)
    solution = Vector{Int32}(undef, n)
    selected_heuristic = Vector{UInt8}(undef, selected_capacity)
    history_objective_values = Vector{Cdouble}(undef, history_capacity)
    history_times_seconds = Vector{Cdouble}(undef, history_capacity)
    heuristic_string = _mqlib_cstring(heuristic)
    hhdata_dir = _mqlib_cstring(_mqlib_hyperheuristic_data_dir())

    input = _MQLibCQUBOInput(
        _MQLIB_C_ABI_VERSION,
        Int32(n),
        pointer(linear),
        Int64(length(quadratic_value)),
        pointer(quadratic_i),
        pointer(quadratic_j),
        pointer(quadratic_value),
        _MQLIB_C_INDEX_BASE_ONE,
        pointer(heuristic_string),
        Cdouble(run_time_limit),
        Int32(random_seed),
        pointer(hhdata_dir),
    )
    result = _MQLibCQUBOResult(
        _MQLIB_C_ABI_VERSION,
        0.0,
        0.0,
        pointer(solution),
        Int32(length(solution)),
        pointer(selected_heuristic),
        Int32(length(selected_heuristic)),
        pointer(history_objective_values),
        pointer(history_times_seconds),
        Int32(length(history_objective_values)),
        0,
    )

    input_ref = Ref(input)
    result_ref = Ref(result)
    status = GC.@preserve linear quadratic_i quadratic_j quadratic_value solution selected_heuristic history_objective_values history_times_seconds heuristic_string hhdata_dir begin
        ccall(
            (:mqlib_solve_qubo, _mqlib_library()),
            Cint,
            (Ref{_MQLibCQUBOInput}, Ref{_MQLibCQUBOResult}),
            input_ref,
            result_ref,
        )
    end

    return (
        status = status,
        result = result_ref[],
        solution,
        selected_heuristic,
        history_objective_values,
        history_times_seconds,
    )
end

function _mqlib_selected_heuristic(buffer::Vector{UInt8}, length::Integer)
    if length <= 1
        return ""
    end

    return String(buffer[1:(Int(length) - 1)])
end

function _print_header(silent::Bool, heuristic::Union{String,Nothing})
    if !silent
        heuristic = something(heuristic, "Hyper-Heuristic")

        print(
            """
            ▷ MQLib
            ▷ Heuristic: $(heuristic)
            ┌────────┬─────────────┬──────────┐
            │  iter  │    value    │   time   │
            ├────────┼─────────────┼──────────┤
            """
        )
    end

    return nothing
end

function _print_footer(silent::Bool)
    if !silent
        println(
            """
            └────────┴─────────────┴──────────┘
            """
        )
    end

    return nothing
end

function _print_iter(silent::Bool, i::Integer, λ::Vector{Float64}, t::Vector{Float64})
    if !silent
        for (λ̄, t̄) in zip(λ, t)
            if isnothing(i)
                @printf("│        │ %11.3f │ %8.2f │\n", λ̄, t̄)
            else
                @printf("│ %6d │ %11.3f │ %8.2f │\n", i, λ̄, t̄)

                i = nothing
            end
        end
    end

    return nothing
end

function _mqlib_args(;
    file_path::String,
    heuristic::Union{String,Nothing},
    random_seed::Union{Integer,Nothing},
    run_time_limit::Float64,
)
    heur = if isnothing(heuristic)
        `-hh`
    else
        `-h $heuristic`
    end

    seed = if isnothing(random_seed)
        ``
    else
        `-s $random_seed`
    end

    return `$heur -fQ $file_path -r $run_time_limit -nv -ps $seed`
end

function unset_heuristic(model)
    set_heuristic(model, nothing)

    return nothing
end

function set_heuristic(model, heuristic::Union{String,Nothing} = nothing)
    MOI.set(model, MQLib.Heuristic(), heuristic)

    return nothing
end

function get_heuristic(model)
    return MOI.get(model, MQLib.Heuristic())
end

function heuristics()
    return sort!(collect(keys(_HEURISTICS)))
end

function show_heuristics()
    for heuristic in heuristics()
        println("$(heuristic):\n  $(_HEURISTICS[heuristic])")
    end

    return nothing
end

end # module
