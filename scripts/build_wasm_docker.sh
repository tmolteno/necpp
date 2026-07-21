#!/bin/bash
# Build the nec2++ WASM target inside the Emscripten Docker image, for hosts
# that don't have a local emsdk install.
#
# Produces nec2pp.js + nec2pp.wasm in the repo root (matching the legacy
# Makefile output names so existing consumers are unaffected).
#
# Replaces the `make wasm` target from the old hand-written Makefile.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="build-wasm"
WASM_IMAGE="emscripten/emsdk:4.0.7"

cd "$PROJECT_DIR"

echo "=== Building WASM via Emscripten Docker image: $WASM_IMAGE ==="
docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "$PROJECT_DIR:/src" \
    -w /src \
    "$WASM_IMAGE" \
    bash -c "
        set -euo pipefail
        emcmake cmake -B $BUILD_DIR -S . \
            -DCMAKE_BUILD_TYPE=Release \
            -DNECPP_BUILD_WASM=ON \
            -DNECPP_BUILD_TESTS=OFF \
            -DBUILD_SHARED_LIBS=OFF
        cmake --build $BUILD_DIR -j\$(nproc)
        # Emscripten emits nec2pp.js + nec2pp.wasm in the build dir; copy them
        # to the repo root for parity with the old Makefile output location.
        cp -f $BUILD_DIR/src/nec2pp.js $BUILD_DIR/src/nec2pp.wasm . 2>/dev/null || true
    "

echo "=== WASM build complete: nec2pp.js + nec2pp.wasm ==="
ls -la nec2pp.js nec2pp.wasm
