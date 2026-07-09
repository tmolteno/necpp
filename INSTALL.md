# Installation from Source

If nec2++ isn't available precompiled in your system, here are the instructions for
building from source code.

## Quick Build (no autotools needed)

Eigen 3.4.0 is bundled in `src/eigen3/` — no external libraries required. Just a
C++17 compiler and `make`:

    sudo apt install g++ make      # Debian/Ubuntu
    make -j4
    sudo make install

Options:

    make DEBUG=1        # bounds-checked debug build
    make TYPECHECK=1    # typesafe integer checking

## macOS Notes

If you see a linker error like "file was built for archive which is not the
architecture being linked", this is usually caused by Homebrew's binutils
conflicting with the system toolchain. Run:

    brew unlink binutils

then rebuild. You can re-link binutils afterwards with `brew link binutils`.

### Building with Xcode

To use necpp in an Xcode project:

1. Add all `.cpp` and `.h` files from `src/` to your project (except
   `nec2cpp.cpp` which contains its own `main()`).
2. Add `src/eigen3/` to your header search paths (`-I src/eigen3`).
3. Create a minimal `config.h`:
   ```c
   #define VERSION "2.1.1"
   #define BUILD_DATE "unknown"
   ```
4. Set C++ Language Dialect to C++17.
5. Link with `-lm` (libm). No other external libraries are needed.

## Compiling for a specific architecture

For best performance, compile with optimizations tuned to your CPU:

    make CXXFLAGS="-std=c++17 -O3 -march=native -Wall -Wextra -Wshadow"

## About Eigen

Nec2++ 2.0.0+ uses **Eigen 3.4.0** (bundled in `src/eigen3/`) for all linear
algebra — safe arrays, 3-vectors, and LU decomposition. No external Eigen, LAPACK,
or BLAS installation is required. The old `--with-eigen`, `--with-lapack`, and
`--with-atlas` configure options have been removed.

## Building for Windows

### Cross-compiling with MinGW

Cross-compile from Linux using the MinGW toolchain:

    sudo apt install g++-mingw-w64-x86-64      # Debian/Ubuntu
    make CXX=x86_64-w64-mingw32-g++

The executable `nec2++.exe` will be built in the project root.

Older 32-bit MinGW:

    make CXX=i686-w64-mingw32-g++

### Compiling with Visual Studio

NEC2++ has been tested with Microsoft Visual Studio 2013.
A solution and project files are provided in the `win32/nec2++/` directory.
Open `nec2++.sln` in Visual Studio and build.

To use necpp as a library in Visual Studio, add the `.cpp` and `.h` files
from `src/` (except `nec2cpp.cpp`), add `src/eigen3/` to your include paths,
and link with `libm`.

## Notes for embedding in other projects

Some source files are header-only and have no corresponding `.cpp` file.
If your build system expects a `.cpp` for every `.h`, the following files
are intentionally header-only and should be included directly:

* `math_util.h` — math constants and the `nec_3vector` class
* `nec_wire.h` — wire intersection geometry
* `safe_array.h` — bounds-checked array template
* `src/eigen3/` — bundled Eigen 3.4.0 headers (include via `-I src/eigen3`)
