#!/bin/bash
# Build the Eigen-based version (current master)
set -e

SRC_DIR="/home/tim/github/necpp/src"
BUILD_DIR="/home/tim/github/necpp/build"
EIGEN_DIR="$SRC_DIR/eigen3"
OUT_DIR="$BUILD_DIR/eigen_build"

rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

cp "$BUILD_DIR/config_eigen.h" "$BUILD_DIR/config.h"

CXX="g++"
CXXFLAGS="-std=c++17 -O0 -g3 -Wall -Wextra -Wshadow -DNEC_ERROR_CHECK=1"
CXXFLAGS="$CXXFLAGS -Dnec_build_date=\"2025-06-25\""
INCLUDES="-I $SRC_DIR -I $BUILD_DIR -isystem $EIGEN_DIR"

echo "=== Building Eigen version (v2.0.0) ==="

SOURCES="c_evlcom c_geometry c_ggrid c_plot_card electromag libNEC matrix_algebra misc nec_context nec_exception nec_ground nec_output nec_radiation_pattern nec_results nec_structure_currents XGetopt"

for src in $SOURCES; do
    echo "  Compiling $src.cpp..."
    $CXX $CXXFLAGS $INCLUDES -c "$SRC_DIR/$src.cpp" -o "$OUT_DIR/$src.o"
done

# Build nec2++ binary
echo "  Linking nec2++..."
ALL_OBJS=""
for src in $SOURCES; do
    ALL_OBJS="$ALL_OBJS $OUT_DIR/$src.o"
done
$CXX $CXXFLAGS $INCLUDES -o "$OUT_DIR/nec2++" "$SRC_DIR/nec2cpp.cpp" $ALL_OBJS -lm -lstdc++

# Build test runner - use linker flag to allow multiple Catch definitions
echo "  Compiling test harness files..."
TEST_FILES="safe_array_tb matrix_algebra_tb c_geometry_tb c_evlcom_tb nec_context_tb"
for tb in $TEST_FILES; do
    echo "    Compiling $tb.cpp..."
    $CXX $CXXFLAGS $INCLUDES -c "$SRC_DIR/$tb.cpp" -o "$OUT_DIR/$tb.o"
done

# Create a separate main for Catch2
cat > "$OUT_DIR/test_main.cpp" << 'EOF'
#define CATCH_CONFIG_MAIN
#include "catch.hpp"
EOF
$CXX $CXXFLAGS $INCLUDES -c "$OUT_DIR/test_main.cpp" -o "$OUT_DIR/test_main.o"

echo "  Linking test_runner..."
$CXX $CXXFLAGS $INCLUDES -Wl,--allow-multiple-definition -o "$OUT_DIR/test_runner" $ALL_OBJS \
    $OUT_DIR/safe_array_tb.o $OUT_DIR/matrix_algebra_tb.o $OUT_DIR/c_geometry_tb.o \
    $OUT_DIR/c_evlcom_tb.o $OUT_DIR/nec_context_tb.o $OUT_DIR/test_main.o -lm -lstdc++

echo "=== Eigen build complete ==="
ls -la "$OUT_DIR/nec2++" "$OUT_DIR/test_runner"
