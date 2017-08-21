# Using nec2++ from C/C++ code

This directory contains an example showing how to call nec2++ from within C or C++ to simulate
the performance of a wire structure.

### Author

Tim Molteno <tim@physics.otago.ac.nz>

## C example

This example requires that the libnecpp library is installed. This is done, either by installing the package in your distribution,
or by compiling the source code (in the directory above this one)

    make test_c

## C++ example

This requires the LAPACK libraries to be installed. These are installable on Debian derived linux using the package libatlas-base-dev.

    make test_cpp
