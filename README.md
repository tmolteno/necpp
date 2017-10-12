# NEC2++ Numerical Electromagnetic Code in C++

This is a free (GPL) electromagnetic simulation software compatable with NEC-2. It has been rewritten from the ground up.

Nec2++ consists of a library that can be called from C++, C, python and Ruby, and so it can incorporated 
into other projects like GUI tools and automatic antenna optimization systems.

There is also an executable necpp that can read antenna description files (like the original). 
Nec2++ is developed on Debian linux, but will work on a variety of other operating systems.

## Features

* NEC-2 syntax compatable.
* C, C++, Python and Ruby libraries included. Ideal for antenna optimization.
* Large designs can be simulated with tens of thousands of elements (to the limit of 64-bit address space)
* Geometry error detection. Throws exceptions if wires intersect or lie too close to one another.
* Simulate in different media (for example antennas in seawater) by modifying the dielectric properties.
* Uses fast numerical routines (BLAS and LAPACK). Can use the Intel MKL or OpenBLAS.

## Citing NEC2++

If you use nec2++, please cite it as follows:

Timothy C.A. Molteno, ''NEC2++: An NEC-2 compatible Numerical Electromagnetics Code'', Electronics Technical Reports No. 2014-3, ISSN 1172-496X, October
2014.

## Documentation

Online documentation built form the source code is available at http://tmolteno.github.io/necpp/. 
A guide to [using nec2++ from python](http://astroelec.blogspot.co.nz/2015/05/modeling-antennas-in-python-with-nec2.html).

## Installation

Instructions are in the file INSTALL.md.

## Links

* http://tmolteno.github.io/necpp/ Documentation
* https://github.com/lncgomz/UCNEC Java-Based GUI project
* https://github.com/tmolteno/python-necpp/ Python packages
* https://github.com/rcnlee/necpp.jl Julia wrapper for the python packages
