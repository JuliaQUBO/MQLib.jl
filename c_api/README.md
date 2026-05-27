# MQLib C ABI

This directory defines the native C ABI for solving QUBO instances through
MQLib without invoking the command-line executable.

The files here are intended to be compiled with the upstream MQLib C++ sources
and packaged as a shared library by `MQLib_jll`. The Julia wrapper does not call
this ABI yet; that integration belongs to the follow-up shared-library issue.

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

## Build Sketch

The exact build recipe belongs in `MQLib_jll`, but a local syntax check against
an upstream MQLib checkout looks like:

```sh
c++ -std=c++11 \
  -Ic_api/include \
  -I/path/to/MQLib/include \
  -fsyntax-only \
  c_api/src/mqlib_c_api.cpp
```

When building the shared library, compile `c_api/src/mqlib_c_api.cpp` together
with the upstream MQLib implementation sources needed by the heuristics. The
upstream executable entry point (`src/main.cpp`) should not be part of the
shared-library target.

The small native example in `c_api/examples/solve_qubo.c` demonstrates both an
explicit QUBO heuristic and the hyperheuristic call shape.
