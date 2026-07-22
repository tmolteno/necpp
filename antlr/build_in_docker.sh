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

# Generate the build config header that common.h / nec2cpp.cpp require.
# The canonical NECPP_VERSION lives in the top-level CMakeLists.txt
# (project(necpp VERSION ...)); read it from there so the two builds cannot
# drift apart. Written to ../build/config.h so the -I ../build flag used in
# the g++ invocations below picks it up. The project() declaration spans two
# lines, so tr flattens the file to one line before sed extracts the version.
NECPP_VERSION=$(tr '\n' ' ' < ../CMakeLists.txt \
  | sed -n 's/.*project(necpp[[:space:]]\{1,\}VERSION[[:space:]]\{1,\}\([0-9][0-9.]*\).*/\1/p' \
  | head -1)
NECPP_BUILD_DATE=$(date +"%Y-%m-%d")
mkdir -p ../build
{
  echo '#ifndef CONFIG_H'
  echo '#define CONFIG_H'
  echo "#define NECPP_VERSION \"${NECPP_VERSION}\""
  echo "#define NECPP_BUILD_DATE \"${NECPP_BUILD_DATE}\""
  echo '#endif'
} > ../build/config.h
echo "  config.h: NECPP_VERSION=${NECPP_VERSION} NECPP_BUILD_DATE=${NECPP_BUILD_DATE}"

echo "=== Compiling full NEC parser ==="
# Compile nec2cpp.cpp separately with main() renamed
mkdir -p build
docker run $DOCKER_OPTS "$IMAGE" \
  g++ -std=c++17 ${CXXFLAGS:--O0 -g} \
    -I . -I ../src -isystem ../src/eigen -I ../build \
    -Dmain=nec2cpp_main_renamed -c ../src/nec2cpp.cpp -o build/nec2cpp.o

docker run $DOCKER_OPTS "$IMAGE" \
  g++ -std=c++17 ${CXXFLAGS:--O0 -g} \
    -I generated -I . -I ../src -isystem ../src/eigen -I ../build \
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
    -l:libantlr4-runtime.a -lm -lstdc++ \
    -o nec_parse

echo "=== Build complete: ./nec_parse ==="
