# Testing nec2++

First create a debug build using the build_debug.sh script

    sh build_debug.sh

You may need to make distclean in the parent directory before this is done.

Then

    make -f Makefile.test

    ## Running a specific test

To run a particular test file

    make DO_TESTS=data/plane_wave_excitation.nec
