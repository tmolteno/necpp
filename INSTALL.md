# Installation from Source

If nec2++ isn't available precompiled in your system, here are the instructions for
building from source code.

## Pre-requisites

Nec2++ requires a C++ compiler and uses the GNU autoconf packages for
keeping track of dependencies. On Debian or a derivative, you can
install these with the following command.

    aptitude install g++ make automake autoconf libtool libatlas-base-dev

## Installation Steps

  1. Install the autoconf and libtool packages.
     On Debian this is done with 
       `aptitude install automake autoconf libtool'
     on other systems you will have to find the appropriate way to do this.

  2. Cenerate the ./configure script. To do this, type
       `make -f Makefile.git'

  3. Then do the usual thing 
       ./configure --without-lapack
       make -j 4
       sudo make install

  4. To use LAPACK, you should install the appropriate atlas system
     aptitude install libatlas-base-dev. And then to the usual thing, i.e.,
        ./configure 
        make -j 4
        sudo make install

## Compiling for a specific architecture

Particularly using gcc-4.0 it is important to specify the architecture you
are using (this will produce approximately a 50% speedup on the athlon-xp
for example. The following configure options will do this:

  1. ./configure CXX=g++-4.0 CXXFLAGS="-O3 -march=athlon-xp"
  2. make
  3. make install

## [Testing only] Using Eigen

Eigen is a fast matrix library for C++. Nec2++ can be build to use eigen rather than its own safe_array class for vectors.

    sudo aptitude install libeigen3-dev
    ./configure --with-eigen

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

