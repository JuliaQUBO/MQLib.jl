module MQLib

import MQLib_jll
using Anneal
using Printf

export get_heuristic, set_heuristic, unset_heuristic, list_heuristics

const _MQLIB_HEURISTICS = Dict{String,String}()

function __init__()
    MQLib_jll.MQLib() do exe
        ms = eachmatch(r"([a-zA-Z0-9]+)\r?\n  ([^\r\n]+)\r?\n?", read(`$exe -l`, String))

        for m in ms
            push!(_MQLIB_HEURISTICS, m[1] => m[2])
        end
    end
end

Anneal.@anew Optimizer begin
    name = "MQLib"
    sense = :max
    domain = :bool
    attributes = begin
        RandomSeed["seed"]::Union{Integer,Nothing} = nothing
        NumberOfReads["num_reads"]::Integer = 1
        Heuristic["heuristic"]::Union{String,Nothing} = nothing
    end
end

function Anneal.sample(sampler::Optimizer{T}) where {T}
    α = QUBOTools.scale(sampler)
    β = QUBOTools.offset(sampler)

    num_reads      = MOI.get(sampler, MQLib.NumberOfReads())
    silent         = MOI.get(sampler, MOI.Silent())
    heuristic      = MOI.get(sampler, MQLib.Heuristic())
    random_seed    = MOI.get(sampler, MQLib.RandomSeed())
    time_limit_sec = MOI.get(sampler, MOI.TimeLimitSec())

    if num_reads <= 0
        error("Number of reads must be a positive integer")
    end

    if !isnothing(heuristic) && !haskey(_MQLIB_HEURISTICS, heuristic)
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

    samples  = Anneal.Sample{T,Int}[]
    metadata = Dict{String,Any}("time" => Dict{String,Any}())

    mktempdir() do path
        qubo_file = joinpath(path, "file.qubo")

        fmt = QUBOTools.QUBO(QUBOTools.BoolDomain(), QUBOTools.MQLibStyle())

        args = _mqlib_args(;
            qubo_file      = qubo_file,
            heuristic      = heuristic,
            random_seed    = random_seed,
            run_time_limit = run_time_limit,
        )

        QUBOTools.write_model(qubo_file, sampler, fmt)

        MQLib_jll.MQLib() do exe
            cmd = `$exe $args`
            
            _print_header(silent, heuristic)
            
            t = 0.0

            for i = 1:num_reads
                lines = readlines(cmd)
                info  = split(lines[begin], ',')

                λ = parse(T, info[4])
                ψ = parse.(Int, split(lines[end], ' '))
                s = Sample{T}(ψ, α * (λ + β))

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

    return Anneal.SampleSet{T}(samples, metadata)
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
    qubo_file::String,
    heuristic::Union{String,Nothing},
    random_seed::Union{Integer,Nothing},
    run_time_limit::Float64,
)
    args = `-fQ $qubo_file -r $run_time_limit -nv -ps`

    if !isnothing(random_seed)
        args = `$args -s $random_seed`
    end

    if isnothing(heuristic)
        args = `$args -hh`
    else
        args = `$args -h $heuristic`
    end

    return args
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
    return sort!(collect(keys(_MQLIB_HEURISTICS)))
end

function show_heuristics()
    for heuristic in heuristics()
        println("$(heuristic): \n  $(_MQLIB_HEURISTICS[heuristic])")
    end

    return nothing
end

end # module