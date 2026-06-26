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

## Autotools Build (traditional)

### Pre-requisites

On Debian or a derivative, you can install the build tools with:

    aptitude install g++ make automake autoconf libtool

### Installation Steps

  1. Install the autoconf and libtool packages.
     On Debian:  `aptitude install automake autoconf libtool`

  2. Generate the `./configure` script:
        `make -f Makefile.git`

  3. Configure, build, and install:
        `./configure`
        `make -j 4`
        `sudo make install`

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
   #define VERSION "2.0.0"
   #define BUILD_DATE "unknown"
   ```
4. Set C++ Language Dialect to C++17.
5. Link with `-lm` (libm). No other external libraries are needed.

## Compiling for a specific architecture

For best performance, compile with optimizations tuned to your CPU:

    make CXXFLAGS="-std=c++17 -O3 -march=native -Wall -Wextra -Wshadow"

With autotools:

    ./configure CXXFLAGS="-O3 -march=native"
    make

## About Eigen

Nec2++ 2.0.0+ uses **Eigen 3.4.0** (bundled in `src/eigen3/`) for all linear
algebra — safe arrays, 3-vectors, and LU decomposition. No external Eigen, LAPACK,
or BLAS installation is required. The old `--with-eigen`, `--with-lapack`, and
`--with-atlas` configure options have been removed.

## Building for Windows

The MinGW toolset (a free compiler for Windows based on GCC) can be used to
compile nec2++ for windows operating systems. This is easily done from a
command line (cygwin shell is best).

  1. ./configure --host=i586-mingw32msvc
  2. make
  3. The executable nec2++.exe can now be found in the src subdirectory.

### Compiling with Visual Studio 2013

NEC2++ has been tested with Microsoft Visual Studio 2013.
Step-by-step instructions

* Build the project inside the win32 subdirectory with Visual Studio 2103.

## Notes for embedding in other projects

Some source files are header-only and have no corresponding `.cpp` file.
If your build system expects a `.cpp` for every `.h`, the following files
are intentionally header-only and should be included directly:

* `math_util.h` — math constants and the `nec_3vector` class
* `nec_wire.h` — wire intersection geometry
* `safe_array.h` — bounds-checked array template
* `src/eigen3/` — bundled Eigen 3.4.0 headers (include via `-I src/eigen3`)
