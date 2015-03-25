# NEC2++ Numerical Electromagnetic Code in C++

This is a free (GPL) electromagnetic simulation software compatable with NEC-2. It has been rewritten from the ground up.

Nec2++ consists of a library that can be called from C++, C, python and Ruby, and so it can incorporated 
into other projects like GUI tools and automatic antenna optimization systems.

There is also an executable necpp that can read antenna description files (like the original). 
Nec2++ is developed on Debian linux, but will work on a variety of other operating systems.

## Features

* Nec2 compatable
* C, C++, Python and Ruby libraries included. Ideal for antenna optimization.
* Large designs can be simulated with tens of thousands of elements (to the limit of 64-bit address space)
* Geometry error detection. Throws exceptions if wires intersect or lie too close to one another.
* Simulate in different media (for example antennas in seawater) by modifying the dielectric properties.
* Uses fast numerical routines (BLAS and LAPACK). Can use the Intel MKL or OpenBLAS.

## Citing NEC2++

If you use nec2++, please cite it as follows:

Timothy C.A. Molteno, ''NEC2++: An NEC-2 compatible Numerical Electromagnetics Code'', Electronics Technical Reports No. 2014-3, ISSN 1172-496X, October
2014.


## Instructions for Linux

nec2++ is available precompiled as part of most modern linux distributions, including Debian, Ubuntu and Fedora.

    sudo aptitude install necpp

## Compiling on Linux

!NOTE! nec2++ now requires LAPACK (unless you use the  --without-lapack configure option)

See the INSTALL file. But here is the short version

    sudo aptitude install libatlas-base-dev
    make -f Makefile.git
    ./configure
    make
    sudo make install

To build a debugging version use

    ./configure --with-bounds

## Instructions for Compiling on Windows

Versions of nec2++ since 1.2.3 now compile fine with the MinGW (http://www.mingw.org/) free compiler tools. 
Just download the source distribution, and follow the unix installation guide (./configure and make).

### Compiling with Visual Studio 2013

NEC2++ has been tested with Microsoft Visual Studio 2013.
Step-by-step instructions

* Build the project inside the win32 subdirectory with Visual Studio 2103.
