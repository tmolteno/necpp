#####################################################################3
#
#	FORTRAN makefile for nec2
#
#	This makefile builds a slightly modified version of the original
#	FORTRAN version of nec2. This is used for the testharness so ensure
#	that nec2++ remains accurate as we develop it.
#
#	The modification is to allow input from stdin and output to stdout.
#	these modifications might be useful to others too.
#
#	Requires g77.
#
#	Author: T Molteno.
#
#
#

#
#	Note about optimization flags and warnings. The compiler issues many many warnings
#	about variables that might be used uninitialized. These may be correct. I have not
#	fixed them since my goal is to simply compare the results of nec2++ with the original
#	FORTRAN. Some of these might however be causing real problems, so perhaps it would be
#	better to fix them.
#

FORTRAN=gfortran

nec2:	nec2dx.f
	${FORTRAN} -Wall -O1 -o nec2 nec2dx.f

somnec:	somnec2d.f
	${FORTRAN} -Wall -O1 -o somnec somnec2d.f

clean:
	rm nec2
	rm somnec

test:
	./nec2 ../test_data/example2.nec ../test_data/example2.outf
