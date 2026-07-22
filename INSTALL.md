# Installation from Source

If nec2++ isn't available precompiled in your system, here are the instructions for
building from source code.

## Quick Build

nec2++ uses [CMake](https://cmake.org) (≥ 3.16). Eigen 5.0.1 is bundled in
`src/eigen/` — no external libraries are required, just a C++17 compiler and
CMake:

    cmake -B build
    cmake --build build -j4
    sudo cmake --install build

This builds the `nec2++` and `nec2diff` command-line tools, the `libnecpp`
static and shared libraries, and installs headers, man pages, and a
`necpp.pc` pkg-config file.

### Build options

Pass these on the configure step (e.g. `-DNECPP_BUILD_TESTS=OFF`):

| Option | Default | Effect |
|--------|---------|--------|
| `CMAKE_BUILD_TYPE` | `Release` | `Debug` adds `-O0 -g`; use `Release` for `-O2`. |
| `NECPP_BUILD_TESTS` | `ON` | Build the Catch2 unit tests (run with `ctest --test-dir build`). |
| `NECPP_BUILD_WASM` | `OFF` | Build the Emscripten/WASM target (needs `emcmake`). |
| `BUILD_SHARED_LIBS` | `ON` | Also build `libnecpp.so` in addition to `libnecpp.a`. |

### Building for debug

    cmake -B build -DCMAKE_BUILD_TYPE=Debug
    cmake --build build -j4

Bounds checking is enabled in the debug build via `NEC_ERROR_CHECK=1` in the
test target.

## Running the tests

    cmake --build build          # builds the test runner
    ctest --test-dir build --output-on-failure

## Cross-compiling

### macOS

CMake works out of the box with Apple Clang:

    cmake -B build
    cmake --build build -j4

To generate an Xcode project instead:

    cmake -B build -G Xcode

If you see a linker error like "file was built for archive which is not the
architecture being linked", this is usually caused by Homebrew's binutils
conflicting with the system toolchain. Run `brew unlink binutils` then rebuild.
You can re-link binutils afterwards with `brew link binutils`.

### Windows

#### Cross-compiling from Linux with MinGW

    sudo apt install g++-mingw-w64-x86-64      # Debian/Ubuntu
    cmake -B build -DCMAKE_C_COMPILER=x86_64-w64-mingw32-gcc \
                   -DCMAKE_CXX_COMPILER=x86_64-w64-mingw32-g++
    cmake --build build

Older 32-bit MinGW:

    cmake -B build -DCMAKE_CXX_COMPILER=i686-w64-mingw32-g++

#### Visual Studio

CMake's Visual Studio generator replaces the old hand-maintained `.vcxproj`
files:

    cmake -B build -G "Visual Studio 17 2022"
    cmake --build build --config Release

### WASM (Emscripten)

Two ways:

**Native emsdk** (if you have it installed):

    emcmake cmake -B build-wasm -DNECPP_BUILD_WASM=ON -DNECPP_BUILD_TESTS=OFF
    cmake --build build-wasm -j4

**Docker wrapper** (no local emsdk needed):

    ./scripts/build_wasm_docker.sh

Both produce `nec2pp.js` + `nec2pp.wasm` exposing a C API
(`nec_create_context`, `nec_process_input`, `nec_get_output`, …).

## Using the library from another project

After `cmake --install`, necpp exposes both pkg-config and a CMake package
config. Pick whichever your build system uses.

**pkg-config:**

    pkg-config --cflags --libs necpp

**CMake (`find_package`):**

```cmake
find_package(necpp REQUIRED)
target_link_libraries(my_app PRIVATE necpp::necpp)
```

## Packaging (DEB / RPM / tarball)

CMake's CPack produces split packages (runtime library, development files,
command-line tools) replacing the old `debian/` directory:

    cd build
    cpack -G DEB    # libnecpp2, libnecpp-dev, necpp
    cpack -G RPM    # needs rpmbuild
    cpack -G TXZ    # plain tarball

## About Eigen

nec2++ uses **Eigen 5.0.1** (bundled in `src/eigen/`) for all linear algebra —
safe arrays, 3-vectors, and LU decomposition. No external Eigen, LAPACK, or
BLAS installation is required. The old `--with-eigen`, `--with-lapack`, and
`--with-atlas` configure options have been removed.

## The ANTLR parser (`nec_parse`)

The experimental ANTLR 4 NEC parser is built separately (it needs Java + the
antlr4 runtime, which are not part of the core CMake build). See
`antlr/README.md`:

    make -C antlr          # builds nec_parse via the Emscripten Docker image

## Notes for embedding in other projects

Some source files are header-only and have no corresponding `.cpp` file.
If your build system expects a `.cpp` for every `.h`, the following files
are intentionally header-only and should be included directly:

* `math_util.h` — math constants and the `nec_3vector` class
* `nec_wire.h` — wire intersection geometry
* `safe_array.h` — bounds-checked array template
* `src/eigen/` — bundled Eigen 5.0.1 headers (include via `-I src/eigen`)
