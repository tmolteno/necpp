# Using nec2++ from C/C++ code

This directory contains an example showing how to call nec2++ from within C or C++ to simulate
the performance of a wire structure.

### Author

Tim Molteno <tim@physics.otago.ac.nz>

## Original Interface

There is an example file called 'example1.nec' that contains the traditional NEC2 cards. You can run this using nec2++ as follows

    nec2++ -i example1.nec -o example1.out

## C example

This example requires that the libnecpp library is installed. This is done, either by installing the package in your distribution,
or by compiling the source code (in the directory above this one)

    make test_c

## C++ example

This requires the LAPACK libraries to be installed. These are installable on Debian derived linux using the package libatlas-base-dev.

    make test_cpp
