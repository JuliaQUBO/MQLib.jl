# MQLib.jl
[MQLib](https://github.com/MQLib/MQLib) wrapper for JuMP (ft. Anneal.jl)

## Installation
```julia
julia> import Pkg

julia> Pkg.add("MQLib")
```

# Basic Usage
```julia
julia> using JuMP, MQLib

julia> Q = [
         -1  2  2
          2 -1  2
          2  2 -1
       ]

julia> model = Model(MQLib.Optimizer)

julia> @variable(model, x[1:3], Bin)

julia> @objective(model, Max, x' * Q * x)

julia> optimize!(model)
```
