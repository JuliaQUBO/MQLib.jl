# MQLib.jl
[MQLib](https://github.com/MQLib/MQLib) wrapper for JuMP (ft. Anneal.jl)

## Installation
```julia
julia> import Pkg

julia> Pkg.add("MQLib")
```

## Basic Usage
```julia
using JuMP, MQLib

Q = [
   -1  2  2
    2 -1  2
    2  2 -1
]

model = Model(MQLib.Optimizer)

@variable(model, x[1:3], Bin)
@objective(model, Max, x' * Q * x)

optimize!(model)
```

## Selecting Heuristics

This wrapper allows one to access all 39 QUBO and Max-Cut Heuristics provided by [MQLib](https://github.com/MQLib/MQLib).
Selecting the method to be used can be achieved via JuMP's attribute interface:

```julia
JuMP.set_optimizer_attribute(model, "heuristic", "ALKHAMIS1998")
```

or by calling MQLib helper functions:

```julia
MQLib.set_heuristic(model, "ALKHAMIS1998")
```

To list available heuristics and their descriptions, run:

```julia
MQLib.show_heuristics()
```