# Using nec2++ from C/C++ code

This directory contains an example showing how to call nec2++ from within C or C++ to simulate
the performance of a wire structure.

### Author

Tim Molteno <tim@physics.otago.ac.nz>

## Original Interface

There is an example file called 'example1.nec' that contains the traditional NEC2 cards. You can run this using nec2++ as follows

    ../build/src/nec2++ -i example1.nec -o example1.out

## C example

This example requires that the libnecpp library is installed. Build with:

    g++ -std=c++17 -I../src -I../src/eigen -I../build test_nec.c \
        -L../build/src -lnecpp -o test_nec

## C++ example

    g++ -std=c++17 -I../src -I../src/eigen -I../build test_cpp.cpp \
        -L../build/src -lnecpp -o test_cpp
