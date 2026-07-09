#!/bin/bash
# Generate ANTLR 4 C++ code and build the full NEC parser.
# All work is done inside the Docker container.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
IMAGE="necpp-antlr4"

cd "$SCRIPT_DIR"

echo "=== Building Docker image: $IMAGE ==="
docker build -t "$IMAGE" .

DOCKER_OPTS="--rm --user $(id -u):$(id -g) -v $PROJECT_DIR:/workspace -w /workspace/antlr"

echo "=== Generating ANTLR 4 parser (full NEC) ==="
docker run $DOCKER_OPTS "$IMAGE" \
  antlr4 -Dlanguage=Cpp -o generated -visitor NECFull.g4

echo "=== Compiling full NEC parser ==="
# Compile nec2cpp.cpp separately with main() renamed
mkdir -p build
docker run $DOCKER_OPTS "$IMAGE" \
  g++ -std=c++17 -O0 -g \
    -I . -I ../src -isystem ../src/eigen -I ../build/simple \
    -Dmain=nec2cpp_main_renamed -c ../src/nec2cpp.cpp -o build/nec2cpp.o

docker run $DOCKER_OPTS "$IMAGE" \
  g++ -std=c++17 -O0 -g \
    -I generated -I . -I ../src -isystem ../src/eigen -I ../build/simple \
    -I /usr/include/antlr4-runtime \
    generated/NECFullLexer.cpp generated/NECFullParser.cpp \
    generated/NECFullBaseVisitor.cpp \
    nec_test.cpp \
    ../src/c_geometry.cpp ../src/c_ggrid.cpp ../src/libNEC.cpp \
    ../src/nec_exception.cpp ../src/nec_context.cpp ../src/nec_ground.cpp \
    ../src/nec_output.cpp ../src/nec_radiation_pattern.cpp \
    ../src/nec_results.cpp ../src/matrix_algebra.cpp ../src/misc.cpp \
    ../src/electromag.cpp ../src/c_evlcom.cpp ../src/c_plot_card.cpp \
    ../src/nec_structure_currents.cpp ../src/XGetopt.cpp \
    build/nec2cpp.o \
    -lantlr4-runtime -lm -lstdc++ \
    -o nec_parse

echo "=== Build complete: ./nec_parse ==="
