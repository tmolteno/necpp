# Testing nec2++

## Quick start

Build a debug binary, then run the regression harness against the reference
engines:

```bash
# 1. Build nec2++ and the nec2c reference engine
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j$(nproc)

# 2. Build nec2c (C reference)
make -C c_src

# 3. Symlink the binaries where the Makefile expects them
ln -sf ../build/src/nec2++ ../src/nec2++
ln -sf ../build/src/nec2diff ../src/nec2diff

# 4. Run all tests (default: herzian_dipole.nec)
make
```

## Running a specific test

```bash
make DO_TESTS=data/gh_flat_spiral.nec
```

## Adding a test

Drop a `.nec` input file into `data/`. The harness will run it through
`nec2++`, `nec2c`, and FORTRAN (if available), then diff the outputs with
`nec2diff`.

## Available test models

The `data/` directory contains over 45 models covering dipoles, patches,
helices, Sommerfeld ground, arrays, and more.
