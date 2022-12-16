module MQLib

using Anneal
import MQLib_jll

const MQLIB_QUBO_HEURISTICS = Set{String}([
    "ALKHAMIS1998",
    "BEASLEY1998SA",
    "BEASLEY1998TS",
    "GLOVER1998a",
    "GLOVER2010",
    "HASAN2000GA",
    "HASAN2000TS",
    "KATAYAMA2000",
    "KATAYAMA2001",
    "LODI1999",
    "LU2010",
    "MERZ1999CROSS",
    "MERZ1999GLS",
    "MERZ1999MUTATE",
    "MERZ2002GREEDY",
    "MERZ2002GREEDYKOPT",
    "MERZ2002KOPT",
    "MERZ2002ONEOPT",
    "MERZ2004",
    "PALUBECKIS2004bMST1",
    "PALUBECKIS2004bMST2",
    "PALUBECKIS2004bMST3",
    "PALUBECKIS2004bMST4",
    "PALUBECKIS2004bMST5",
    "PALUBECKIS2004bSTS",
    "PALUBECKIS2006",
    "PARDALOS2008",
])

Anneal.@anew Optimizer begin
    name   = "MQLib"
    sense  = :max
    domain = :bool
    attributes = begin
        RandomSeed["seed"]::Union{Integer,Nothing} = nothing
        NumberOfReads["num_reads"]::Integer        = 1
        QUBOHeuristic["heuristic"]::String         = "ALKHAMIS1998"
    end
end

@doc raw"""
    MQLib.Optimizer{T} where T

Available QUBO Heuristics are:

- ALKHAMIS1998: Simulated annealing
- BEASLEY1998SA: Simulated annealing
- BEASLEY1998TS: Tabu search
- GLOVER1998a: Tabu search
- GLOVER2010: Tabu search with long-term memory
- HASAN2000GA: Genetic algorithm
- HASAN2000TS: Tabu search
- KATAYAMA2000: Genetic algorithm with k-opt local search
- KATAYAMA2001: Simulated annealing
- LODI1999: Genetic algorithm
- LU2010: Genetic algorithm with tabu search
- MERZ1999CROSS: Genetic algorithm, with crossover only
- MERZ1999GLS: Genetic algorithm, with crossover and local search
- MERZ1999MUTATE: Genetic algorithm, with mutation only
- MERZ2002GREEDY: GRASP without local search
- MERZ2002GREEDYKOPT: k-opt local search with GRASP
- MERZ2002KOPT: k-opt local search with random restarts
- MERZ2002ONEOPT: 1-opt local search with random restarts
- MERZ2004: Genetic algorithm with k-opt local search
- PALUBECKIS2004bMST1: Tabu search procedure
- PALUBECKIS2004bMST2: Iterated tabu search
- PALUBECKIS2004bMST3: Tabu search with GRASP
- PALUBECKIS2004bMST4: Tabu search with long-term memory
- PALUBECKIS2004bMST5: Iterated tabu search
- PALUBECKIS2004bSTS: Tabu search procedure
- PALUBECKIS2006: Iterated tabu search
- PARDALOS2008: Global equilibrium search
""" Optimizer

function Anneal.sample(sampler::Optimizer{T}) where {T}
    α = QUBOTools.scale(sampler)
    β = QUBOTools.offset(sampler)

    num_reads = MOI.get(sampler, MQLib.NumberOfReads())

    if num_reads <= 0
        error("Number of reads must be a positive integer")
    end

    qubo_heuristic = MOI.get(sampler, MQLib.QUBOHeuristic())
    
    if qubo_heuristic ∉ MQLIB_QUBO_HEURISTICS
        error("Invalid QUBO Heuristic code '$qubo_heuristic'")
    end

    random_seed = MOI.get(sampler, MQLib.RandomSeed())

    if isnothing(random_seed)
        random_seed = trunc(Int, time()) % 65_536
    else
        random_seed = random_seed % 65_536
    end

    time_limit_sec = MOI.get(sampler, MOI.TimeLimitSec())

    run_time_limit = if isnothing(time_limit_sec)
        1.0 / num_reads
    else
        time_limit_sec / num_reads
    end

    runtime = 0.0
    samples = Anneal.Sample{T,Int}[]

    mktempdir() do path
        qubo_file = joinpath(path, "file.qubo")

        write(qubo_file, sampler, QUBOTools.QUBO(; style = :mqlib))

        MQLib_jll.MQLib() do exe
            cmd = `
                $exe
                -h  $qubo_heuristic
                -fQ $qubo_file
                -r  $run_time_limit
                -nv
                -ps
            `

            for _ = 1:num_reads
                let lines = readlines(cmd)
                    info = split(lines[begin], ',')

                    λ = parse(T, info[4])
                    ψ = parse.(Int, split(lines[end], ' '))
                    
                    push!(samples, Sample{T}(ψ, α * (λ + β)))

                    runtime += parse(Float64, info[5])
                end
            end
        end
    end

    metadata = Dict{String,Any}(
        "time" => Dict{String,Any}(
            "effective" => runtime,
        )
    )

    return Anneal.SampleSet{T}(samples, metadata)
end

end # module