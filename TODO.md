# TODO list for nec2++

* Add remaining commands to the API for C, Ruby and Python.
* Test and make sure that the Somerfeld Ground code is working (compare outputs with the FORTRAN version)
* Get CurrentInput comparison going with nec2diff

* Work on Doxygen Source Documentation
* Get rid of libnec's dependance on the C++ standard library. See Below.
* Export only the C API functions from libnec. This requires some extra know-how that I don't have right now -- probably an exports file. 
  When done, this will mean that libnec can be linked to without additional libc++ linking.

## NEC Syntax Improvements

* Handle commenting out lines (Look into how this should be done). Some folk seem to use a ' at the start of the line.
* Add a geometry method (and nec card) for catenary hung wires (cf NEC4) CW card input.


## NEC improvements

* Export radiation patterns (command line switch) as CSV files. This will allow them to be plotted using other plotting tools.
* Add Leeson Corrections for tapered wires (these are frequency dependant).
* Add Monte Carlo simulations that measure the sensitivity of a design to random fluctuations in the geometry. Suggestion from Mike Pot.
* Add output option for the impedance (admittance) matrix to go to the output file.

## Using Eigen

* To avoid dependency on Lapack, we are going to shift the LU decomposition to using the Eigen library.

