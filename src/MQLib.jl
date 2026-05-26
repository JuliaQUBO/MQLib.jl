module MQLib

using Printf

import MQLib_jll
import QUBOTools
import QUBODrivers
import MathOptInterface as MOI

const __VERSION__ = v"0.1.0"
const _HEURISTICS = Dict{String,String}()
const _MQLIB_STDIN = "/dev/stdin"

abstract type _MQLibInput end

struct _MQLibFileInput <: _MQLibInput
    file_path::String
end

struct _MQLibStdinInput <: _MQLibInput
    data::Vector{UInt8}
    file_path::String
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

    V = Set{Int}(1:n)

    model = QUBOTools.Model{Int,T,Int}(
        V, L, Q;
        scale  = α,
        offset = β,
        sense  = :max,
        domain = :bool,
    )

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

    samples  = QUBODrivers.Sample{T,Int}[]
    metadata = Dict{String,Any}(
        "time"   => Dict{String,Any}(),
        "origin" => Dict{String,Any}(
            "name"      => "MQLib",
            "version"   => __VERSION__,
            "heuristic" => heuristic,
        ),
    )

    mktempdir() do temp_path
        input = _mqlib_input(model, joinpath(temp_path, "model.qubo"))

        args = _mqlib_args(;
            file_path = _mqlib_file_path(input),
            heuristic,
            random_seed,
            run_time_limit,
        )

        let exe = MQLib_jll.MQLib()
            cmd = `$exe $args`

            _print_header(silent, heuristic)

            t = 0.0
            
            for i = 1:num_reads
                lines = _mqlib_readlines(cmd, input)
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

            metadata["time"]["effective"] = t
        end
    end

    return QUBOTools.SampleSet{T}(samples, metadata; sense = :max, domain = :bool)
end

function _mqlib_input(model::QUBOTools.AbstractModel, file_path::String)
    if _supports_mqlib_stdin()
        io = IOBuffer()

        QUBOTools.write_model(io, model, QUBOTools.Format{:qubo}(; style = :mqlib))

        return _MQLibStdinInput(take!(io), _MQLIB_STDIN)
    else
        QUBOTools.write_model(file_path, model, QUBOTools.Format{:qubo}(; style = :mqlib))

        return _MQLibFileInput(file_path)
    end
end

_mqlib_file_path(input::_MQLibInput) = input.file_path

_supports_mqlib_stdin() = Sys.isunix() && ispath(_MQLIB_STDIN)

_mqlib_readlines(cmd::Cmd, input::_MQLibFileInput) = readlines(cmd)

function _mqlib_readlines(cmd::Cmd, input::_MQLibStdinInput)
    return readlines(pipeline(cmd; stdin = IOBuffer(input.data)))
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
