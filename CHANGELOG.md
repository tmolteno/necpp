## Unreleased

* Fix error in Sommerfeld ground interpolation thanks to Alexander Schewelew
* Update checking for blas, so that other blas can be used like openblas.

## Vresion 1.7.4
* Fix bug in reading comments longer than 80 characters (issue #47)
* Fixed bug in nec_context::fblock when setting up symmetry_array matrix for plane symmetry.
* Many changes to the 64-bit integer support for large matrices.

## Version 1.7.3
* Changes to allow compilation on CLANG. Change to use std::abs() and friends.
* Minor bug fixes in example c code.
* Minor bug in fallthrough on printing complex numbers.

## Version 1.7.2
* Fixed bug in helix introduced in 1.7.1 - had not updated all calls to the helix function (thanks to Yoshi Takeyasu for pointing this out)

## Version 1.7.1
* Changed API of c_geometry::helix() to put the tag id and segment count at the start.

## Version 1.7.0
* Added new functions for nec_get_gain() to get the radiation pattern elements.
* added nec_medium_parameters to the C-style API. The user can now alter the permittivity and permeability of the environment via the API.
* Added accessor methods to the radiation pattern to get all the theta and phi angles
* Added 2D matrix for radiation pattern gain_. Much easier to plot radiation patterns.
* nec_geometry_complete() now takes only one parameter. The second parameter was unused.
* New c_geometry->reflect and c_geometry->generate_cylindrical_structure(int itx, int nop) helper functions to generate symmetric structures.
* New c_ground methods for finding out the ground parameters
* Split the ground enviroment computation from the output code. 
* Fix make dependencies to allow parallel make
* Allow building with gcc 5.1

## Version 1.6.2 [April 2015]
* Added new functions for nec_excitation_voltage, nec_excitation_current and nec_excitation_planewave to simplfy C-style API calls
* More documentation
* Published python-necpp as a separate project.

## Version 1.6.1 [March 29 2015]
* Added python library for the c-stye API.
* Fixed bug.  "SP" with multiple "SC" doesn't work #10.
* Yoshi Takeyasu. Change instances of int to int64_t, to allow simulations up to 140 GB.
* Bugfix in c_geometry::cmss
* C-style library now supports error messages. This changes the public API.
* Added missing library function nec_kh_card(). 
* Cleaned up bitrot in example code.
* Added C-style support for patch command. nec_patch.
* Added C-style support for GM card (move command). nec_gm_card.
* Added python based testharness
* Start to use matrix operations for multi-dimensional arrays. safe_matrix class.
* Errors are thrown if zeros are specified for ground conductivity and dielectric constant.
* Added Reflection nec_gx_card to C-style interface.
* Added Reflection gx_card to c_geometry class.
* Catch out of memory errors and create nice error messages
* Fix memory allocation errors in patch simulations.

## Version 1.5.4 [December 2014]
* Changed default compiler options to remove debugging symbols
* Fixed bug printing maximum gain on the command line.
* Remove vector test code from main loop 
* Renamed the c_geometry->n to n_segments

## Version 1.5.3 [October 2014]
* Compatability with gcc4.9
* Compatability with gcc4.9
* Improved ANTLR parser.
* Updated PyNEC directory to use numpy instead of numarray. Fixed wrap/*.cpp files [Thanks to James David]
* Moved to using Makefile.git, and updated configure.in to configure.ac
* Include the necpp headers in both CFLAGS and CPPFLAGS

## Version 1.5.2 [May 2013]
* Changed the autoconf system to allow building with the eigen3 package
* Added a genetic optimzer written in Ruby

## Version 1.5.1 [March 2012]
* Improved error messages on intersection. Now reports tag_id as well as wire number.
* Moved source to github
* Added automatic build date (removed their manual definition in common.h)
* New code directory structure. More logical, all testharness files moved to the testharness directory.

## Version 1.5 - 2010-10-02
* Added a test for segment midpoint intersections. If two wires are joined at acute angles to each other, this
  can lead to inaccurate simulations. If either the first or last segment midpoint of a wire lies inside any other wire
  an exception is thrown with a message.
* Better exception handling in the Ruby binding. Exception's now have their messages passed through.
* Improved preporting of power gain. This means that hopefully (not quite fully tested) the gain
  normalization does not any longer have any effect on the get_maximum_gain functions e.t.c.
* Modified to work with atlas 3.8
* Fixed up the benchmark to avoid wires at a very acute angles (violating the midpoint intersection test)

## Version 1.4 - 2009-01-02
* Added accessor functions for radiation pattern statistics.
* Renamed some methods
* Added Ruby Wrapper to the C library
* Added code to detect wire intersections and throw an exception.
* Renamed some methods in the C library interface. This makes a more logical grouping when 
  getting statistics from radiation patterns.
* Added Ruby Wrapper to the C library libnecpp. This code is in the directory Ruby.

## Version 1.3
* Using the correct value for speed of light throughout (was using 299.8 as default freq, and this caused slight 
  problems. Now using 1 meter as default wavelength (rather than 299.8 as the default frequency)
* Added error message when segment length is less that 0.01 wavelengths long.
* Added accessor functions for RHCP and LHCP receiving gain to nec_radiation_pattern
	nec_float get_power_gain_tot(int theta_index, int phi_index) const
	nec_float get_power_gain_rhcp(int theta_index, int phi_index) const
	nec_float get_power_gain_lhcp(int theta_index, int phi_index) const
* Changes geometry parsing to make ix,iy,iz local variables. Set iy=0 before calling reflect for GR card.

## Version 1.2.9
* Improvements to the LAPACK use. Autoconf now checks for the correct functions.

## Version 1.2.8
* We NOW are using LAPACK. Huge speed bonus. Run configure with the option --with-lapack to compile this in. You will need to install some LAPACK libraries. Some work needed in the autoconf files (configure.in) to make this work correctly on more platforms. On debian use
	aptitude install libatlas-base-dev

## Version 1.2.7
* Some work on using ATLAS for the Gaussian Elimination.
* cleaned up some code that produced warnings in gcc 4.2 (const char coercions)  

## Version 1.2.6
* Added output flag for XML -x
* Cleaned up output file handling in nec_results (a little)

## Version 1.2.5

* Fixed a bug in handling of patches. Thanks to the hard work of Neoklis and Jerry Burke (original author),  this has been sorted out. It is tested with neokliks_bug.nec (in test data) and caused a segmentation fault in all versions, C, C++ and FORTRAN.
* Some code has been added for handling arbitrary media - permittivity and permeability. A new card has been added as well.

## Version 1.2.4

* (Remi Sassolas) Major Progress in the Python Port
* (Remi Sassolas) Fixed crashing bug in structure current printing. This bug would happen when someone absent-minded (or just silly) would ask the currents to be printed using the output format designed for a receiving pattern, but wouldn't use an incident plane wave as the excitation type. Then NEC-2 would crash (segmentation fault)
* (Remi Sassolas) New Results object for handling structure currents.


## Version 1.2.3

* Added new test case (pjw_small.nec) that appears to fail on the radiation pattern result when compared to original FORTRAN. This was caused by an error in the original FORTRAN parser. Decided to do better and use a proper LL(k) parser generator.
* Added ANTLR grammar (see http://www.antlr.org) that generates a parser for nec files. This was motivated by a silent failure the caused the wrong result to be printed. For the moment this development is separate (inside the antlr directory). However by version 1.3 we will be using an ANTLR grammar.
* Moved lots of code from header files into .cpp files. This makes the libnecpp more usable as a C++ shared library.
* Added a new method 'calculate_network_data()' to the nec_context class. This separates out the calculations from the print_network_data() method. This code is still ugly but should be functional.
* Moved development to alioth.debian.org to support the python programming effort.
* src/PowerBudget.h: Allow for blank line following POWER INPUT to cope with FORTRAN output
* src/nec2cpp.cpp: Check number of parameters in the TL card and return an error if not all parameters are supplied. This resolves an ambiguity when blank parameters are present (we can't test for these).
* src/nec_radiation_pattern.h: Added accessor functions to the nec_radiation_pattern results class.
* Added an example file "example/test_nec.cpp" to show how to use a results object from a C++ program.

## Version 1.2.2

* Trap condition where n<=1 in c_geometry::connect_segments pending finding where the bug is. The array icon1 does not have its values initialized under these circumstances and bad things happen (walk all over memory).
* Cleaned up the connect_segments code in c_geometry.h. Moved many variable declarations into local contexts.
* Added L1-norm and Euclidian norm to the nec_3vector class (use it in c_geometry::connect_segments)
* Clean up local variables in nec_context::nhfld()
* Clean up local variables in nec_context::pcint()
* Fixed conversion of int to bool when calling nec_context::efld(). Removed unnecessary ij variable.

## Version 1.2.1

* Make ksymp a private member of nec_ground. Replaced it with a boolean accessor method present() that indicates whether the ground is present. Then removed some ugly code in nec_context::efld().
* Switched on all warnings I can think of (more than -Wall) in the debug build. This found a lot of shadowed variables.
* Removed shadowed variables in c_geometry.
* Switched normalized receiving pattern output to use the results object.
	
	
## Version 1.2

* Cleaned up:
	degree to radian conversion. Speedup.
	Radiation Pattern Code. Easier to read (plenty of work left)
* New Results API allows easy access to simulation results, currently includes the following output
	* Antenna Input parameters
	* Normalized Receiving Pattern
	* Radiation Pattern
	* Structure excitation
* Note: The new radiation pattern output will generate phases of complex numbers that are between -180 and +180. For example -195 will appear as +165 degrees. 
* Enabled the writing of some data to standard output (-s command line option)
* NEW: Data can be written in comma-separated-value format for easy importing into other plotting programmes.
* Fixed bug in printing STRUCTURE EXCITATION data (the segment number was incorrect)
* Testharness fixes
* Faster benchmark code (similar output values) This is because we are using some dynamic benchmarking of clients in a heterogenous multiprocessing systemm.


## Version 1.1.4

* Removed the inclusion of c_geometry.h from nec_context.h. This got rid of some unnecessary dependencies. Also changed nec_context to point to a geometry object (rather than own it).
* Added classes for representation of geometry -- to construct XML geometry descriptions.
* Added enum for excitation type
* New nec_3vector class to simplify vast amounts of code.
* Fixed bug in current excitation (inside etmns). Bug exists also in nec2c. This occurs if the excitation type (excite_type) parameter is EXCITATION_CURRENT.
* Changed the parameters of the rp_card function to explicitly take the XNDA parameters as separate integers. This avoids 0500 being passed in as an octal number!
* Code cleanup in c_geometry::tbf()

## Version 1.1.3

* Slight changes to benchmarking output.
* Notes on compilation

## Version 1.1.2

* New header files for the libnecpp.h.
* Moved to libtool for library generation.
* Minor code cleanup in c_evlcom (variable renaming and commenting)
* Added endl() to version.


## Version 1.1.1

21 January 2005.

* Cleaned up and reworked patch subdivision.
* Fixed a bug in interpolate() caused by not having the precalculated data static. Thanks to Neokolis for pointing this out.
* Cleaned up rom2() and rechecked program flow.


## Version 1.1.0

20 January 2005.

* Got LAPACK going. There is generally a 100% speedup. However the LAPACK routines sometimes
produce different answers to the built in LU decomposition routine? Therefore LAPACK LU decompositions is still not included in the default build.
* Fixed up a re-initialization problem in the temporary geometry files.
* Several Improved testharness code.
* Moved temporary geometry into the geometry class. Set this up in the geometry_complete() method.
* Added #include <unistd.h> to misc.cpp
* Ignore Radiation Pattern Polarization angles in the testharness where the power level is -999 dB. These have no physical meaning and are often different with LAPACK.

## Version 1.0.8

December 2004.

* Clean up of nec_context member variables s,b. Renamed m_s, and m_b
* Added member function to c_geometry for testing whenther we should be using a thin wire approximation or not.
* Fix problem introduced in the Sommerfeld-Norton ground condition, added more testharness stuff.

## Version 1.0.7

December 2004.

* Switched command-line option parsing to XGetopt for cross platform compatibility
* Fixed major bug in Win32 executable that caused a crash on file output. Added Visual Studio Project files
* Changed member 'near' of nec_context to m_near to avoid a conflict with the VC++ reserved word 'near'
* Modified the error message macro to allow compilation on Visual C++. This causes error messages to NOT have a content
if nec2++ is compiled with a non-C99 standard compiler.

## Version 1.0.6
December 2004.

* Fix for non-initialization of voltage sources if no excitations were specified (Thanks to Jerome Plet)
* Cleaned up code for bad loop in nec_context (Thanks to Neoklis Kyriazis for this)
* Removed old Numerical Greens Functions code and variables -- these were not doing anything (ib11 e.t.c.)
* All system exit (stop()) calls have been removed and replaced by exceptions (apart from those inside the main command line programme)


## Version 1.0.5 
December 2004.

* Fix for not clearing temporary geometry correctly.
* Moved more stop(-1) commands to throwing nec_exceptions.
* Improved Doxygen Comments
* Added arc and helix commands to the nec_context object.


## Version 1.0.4 
November 2004.

This version includes significant changes. nec2++ can now be called from a C API, although I am still figuring out how
to remove the requirement for linking to the standard C++ library! Help would be appreciated here. See the examples
directory in the source tree.

This version has also moved to KDevelop and autoconf as the primary development environment. This means that nec2++
can be build using the industry standard './configure', 'make' and 'make install' commands.

Please note that my aim is to NOT require KDevelop to build this project, but rather I am using it here. I might
well have missed something when creating the source distribution. If I have, please let me know.
