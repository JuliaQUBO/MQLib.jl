#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="$(mktemp -d)"
trap 'rm -rf "$BUILD_DIR"' EXIT
DEFAULT_MQLIB_UPSTREAM_REF="585496274af5abb0849d0d47e135496b4688680b"
MQLIB_UPSTREAM_REF="${MQLIB_UPSTREAM_REF:-$DEFAULT_MQLIB_UPSTREAM_REF}"

if [[ -n "${MQLIB_UPSTREAM_DIR:-}" ]]; then
    UPSTREAM_DIR="$MQLIB_UPSTREAM_DIR"
else
    UPSTREAM_DIR="$BUILD_DIR/MQLib-upstream"
    git clone --depth 1 https://github.com/MQLib/MQLib.git "$UPSTREAM_DIR"
    if [[ "$MQLIB_UPSTREAM_REF" != "HEAD" ]]; then
        git -C "$UPSTREAM_DIR" fetch --depth 1 origin "$MQLIB_UPSTREAM_REF"
        git -C "$UPSTREAM_DIR" checkout --detach FETCH_HEAD
    fi
fi

echo "Using upstream MQLib $(git -C "$UPSTREAM_DIR" rev-parse --short HEAD)"

if [[ ! -d "$UPSTREAM_DIR/include" || ! -d "$UPSTREAM_DIR/hhdata" ]]; then
    echo "MQLIB_UPSTREAM_DIR must point to an upstream MQLib checkout with include/ and hhdata/" >&2
    exit 2
fi

mapfile -t UPSTREAM_SOURCES < <(find "$UPSTREAM_DIR/src" -name '*.cpp' ! -name main.cpp | sort)

c++ -std=c++11 -fPIC \
    -I"$ROOT/c_api/include" \
    -I"$UPSTREAM_DIR/include" \
    "${UPSTREAM_SOURCES[@]}" \
    "$ROOT/c_api/src/mqlib_c_api.cpp" \
    -shared \
    -o "$BUILD_DIR/libmqlib_c_api.so" \
    -lm

cc -std=c99 \
    -I"$ROOT/c_api/include" \
    "$ROOT/c_api/examples/solve_qubo.c" \
    -L"$BUILD_DIR" \
    -lmqlib_c_api \
    -Wl,-rpath,"$BUILD_DIR" \
    -o "$BUILD_DIR/solve_qubo"

cc -std=c99 \
    -I"$ROOT/c_api/include" \
    "$ROOT/test/c_api_smoke.c" \
    -L"$BUILD_DIR" \
    -lmqlib_c_api \
    -Wl,-rpath,"$BUILD_DIR" \
    -o "$BUILD_DIR/c_api_smoke"

RUN_DIR="$BUILD_DIR/run-without-hhdata"
mkdir -p "$RUN_DIR"

(
    cd "$RUN_DIR"
    "$BUILD_DIR/solve_qubo" "$UPSTREAM_DIR/hhdata"
    "$BUILD_DIR/c_api_smoke" "$UPSTREAM_DIR/hhdata"
)
