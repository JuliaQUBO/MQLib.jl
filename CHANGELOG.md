# Changelog

## Unreleased

## v0.6.0 - 2026-06-06

- Added a Julia `ccall` solve path that uses `MQLib_jll.libmqlib_c_api`.
- Required `MQLib_jll` 0.1.2 or newer so the shared-library product is available from the published JLL.
- Kept the existing executable-backed solve path as a defensive fallback if the C ABI library is unavailable.
- Used the executable fallback for default hyperheuristic solves when the JLL artifact does not include hyperheuristic model files.
- Updated the BinaryBuilder recipe sketch to install MQLib hyperheuristic model files under `share/mqlib/hhdata`.

## v0.5.0 - 2026-05-27

- Raised the Julia compatibility floor from 1.9 to 1.10.
- Updated QUBODrivers compatibility to 0.4 and QUBOTools compatibility to 0.12.
- Updated CI coverage to test Julia 1.10 and latest stable Julia.
- Added a documented C ABI wrapper for solving QUBO instances through MQLib without invoking the command-line executable.
- Added native C ABI smoke tests and contract checks for the exported C API.
- Added a BinaryBuilder recipe sketch for building `MQLib_jll` with both the existing executable product and a future `libmqlib_c_api` shared library product.
