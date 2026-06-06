# MQLib C ABI

This directory defines the native C ABI for solving QUBO instances through
MQLib without invoking the command-line executable.

The files here are intended to be compiled with the upstream MQLib C++ sources
and packaged as a shared library by `MQLib_jll`. The Julia wrapper uses this
ABI through the `libmqlib_c_api` product exported by `MQLib_jll` v0.1.2 and
newer. The executable-backed path remains as a defensive fallback if the library
product is unavailable, or if a default hyperheuristic solve needs model files
that are not packaged in the JLL.

## API

Include `c_api/include/mqlib_c_api.h` from C, C++, or Julia `ccall` bindings.
The public ABI exposes only C scalars, pointers, and structs:

- `MQLibCQUBOInput` passes the problem dimension, linear coefficients, sparse
  off-diagonal quadratic coordinate arrays, heuristic selection, runtime limit,
  and random seed.
- `MQLibCQUBOResult` passes caller-owned buffers for the objective value,
  binary solution vector, runtime, selected heuristic, and optional incumbent
  history.
- `mqlib_solve_qubo` returns a stable `MQLibCStatus` code. It does not transfer
  ownership of any caller buffer.

The wrapper validates ABI inputs before constructing MQLib C++ objects and
reports those failures through status codes. Unrecoverable exits inside
upstream heuristic internals are still upstream behavior until those internals
are made status-code aware.

Quadratic coordinates may be zero-based or one-based, selected by
`index_base`. Diagonal terms belong in `linear`; quadratic arrays must contain
only off-diagonal entries. MQLib interprets off-diagonal entries as the upper
triangle of the symmetric QUBO matrix, matching the existing MQLib `.qubo` file
format.

Set `heuristic` to a QUBO heuristic code such as `ALKHAMIS1998` to run that
heuristic. Set it to `NULL` or an empty string to run the hyperheuristic path.
For hyperheuristic runs, set `hyperheuristic_data_dir` to the directory
containing MQLib's random-forest model files, usually the upstream `hhdata`
directory. If no model files are found, `mqlib_solve_qubo` returns
`MQLIB_STATUS_HYPERHEURISTIC_DATA_NOT_FOUND`.

## BinaryBuilder Recipe

The `MQLib_jll` build recipe in `jll/build_tarballs.jl` builds the existing
`MQLib` executable product and a shared library product named
`libmqlib_c_api`. The library compiles `c_api/src/mqlib_c_api.cpp` with the
upstream MQLib implementation sources, excluding the upstream executable entry
point (`src/main.cpp`), installs this header as `mqlib_c_api.h`, and installs
the hyperheuristic model files under `share/mqlib/hhdata`.

Downstream Julia code should call the library product exported by `MQLib_jll`,
for example:

```julia
ccall((:mqlib_c_abi_version, libmqlib_c_api), Cint, ())
```

The public ABI is versioned by `MQLIB_C_ABI_VERSION`. Any incompatible change
to the structs, status codes, or exported functions must bump that value and be
released through a new `MQLib_jll` build.

A local syntax check against an upstream MQLib checkout looks like:

```sh
c++ -std=c++11 \
  -Ic_api/include \
  -I/path/to/MQLib/include \
  -fsyntax-only \
  c_api/src/mqlib_c_api.cpp
```

The small native example in `c_api/examples/solve_qubo.c` demonstrates both an
explicit QUBO heuristic and the hyperheuristic call shape:

```sh
./solve_qubo /path/to/MQLib/hhdata
```
