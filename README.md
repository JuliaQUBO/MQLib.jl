# MQLib.jl
[![DOI](https://zenodo.org/badge/568213607.svg)](https://zenodo.org/badge/latestdoi/568213607)
[![QUBODRIVERS](https://img.shields.io/badge/Powered%20by-QUBODrivers.jl-%20%234063d8)](https://github.com/JuliaQUBO/QUBODrivers.jl)

[MQLib](https://github.com/MQLib/MQLib) wrapper for JuMP.

MQLib.jl uses the `libmqlib_c_api` shared-library product provided by
`MQLib_jll` v0.1.2 and newer. The executable-backed path remains as a
defensive fallback if the library product is unavailable, or if a default
hyperheuristic solve needs model files that are not packaged in the JLL.

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

This wrapper allows one to access the QUBO and Max-Cut heuristics provided by [MQLib](https://github.com/MQLib/MQLib).
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

## Reads, Seeds, and Time Limits

MQLib.jl supports the standard `QUBODrivers.RandomSeed()` attribute through the
raw `"seed"` optimizer attribute. `QUBODrivers.FinalNumberOfReads()` controls
the number of emitted samples; when it is unset, it falls back to `"num_reads"`.

`MOI.TimeLimitSec()` is divided evenly across the emitted reads and passed to
MQLib as the per-heuristic runtime limit. If no time limit is set, each read
uses a one-second total default split across the requested reads.
