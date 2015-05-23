###########################################################################
#	
#	MAKEFILE for the NEC++ testharness
#	
#	Author: Tim Molteno.
#
#
#	Copyright (C) 2004-2012  Timothy C.A. Molteno
#	
#	This program is free software; you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation; either version 2 of the License, or
#	(at your option) any later version.
#	
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#	
#	You should have received a copy of the GNU General Public License
#	along with this program; if not, write to the Free Software
#	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#

all:	check_build test python

BUILD=../src
TARGET = ${BUILD}/nec2++
NEC2DIFF = ${BUILD}/nec2diff


NEC2C = c_src/nec2c
NECF = FORTRAN/nec2
#
#	This target makes sure that the C and FORTRAN executables have been built
#	Please note that the testharness assumes that you have built the debug
#	version of nec2++ using ./configure --with-bounds
#
check_build: ${NEC2C} ${TARGET}


############################################################################
#	
#	TESTHARNESS
#
#	This tests the output of the C++ port against the nec2c, C 
#	port as well as the original FORTRAN
#
NEC_TESTS := $(shell find data/ -name "example*.nec")
SOM_TESTS := $(shell find data/ -name "sommer*.nec")
ALL_TESTS := $(shell find data/ -name "*.nec")

DO_TESTS = ${NEC_TESTS} ${SOM_TESTS}
#DO_TESTS = ${ALL_TESTS}

############################################################################
#	
#	Rewrite Rules are used to automate the generation
#	and comparison of output files.

# 	Generate an output file from nec2c
%.nec %.outc:
	- ./${NEC2C} -i $*.nec -o temp.$(*F); mv temp.$(*F) $@
	
# 	Generate an output file from nec2cpp
%.nec %.outcpp:
	- ./$(TARGET) -i $*.nec -o $@

# 	Generate an output file from the FORTRAN version
%.nec %.outf:
	- ./$(NECF) $*.nec $@

#	Compare two outputs (c and cpp)
%.nec %.diff: %.outc %.outcpp
	./${NEC2DIFF} $*.outc $*.outcpp > $@
	./${NEC2DIFF} $*.outf $*.outcpp >> $@

#	Compare two outputs (FORTRAN and cpp)
%.nec %.dif: %.out %.outcpp
	./${NEC2DIFF} $*.out $*.outcpp > $@

test:
	echo "testing ${DO_TESTS}"
	rm -f $(DO_TESTS:.nec=.diff)
	rm -f $(DO_TESTS:.nec=.outcpp)
	rm -f $(DO_TESTS:.nec=.outc)
	rm -f $(DO_TESTS:.nec=.outf)
	make test_aux

test_f: ${NECF} $(DO_TESTS:.nec=.outf)
	echo "Fortran test complete"

# 	Change the following line to
#	test_aux:	${NEC2DIFF} nec2cpp $(NEC_TESTS:.nec=.outcpp) $(NEC_TESTS:.nec=.outc) $(NEC_TESTS:.nec=.diff) $(NEC_TESTS:.nec=.dif)
#	in order to include the FORTRAN code in the test suite.

test_aux: $(DO_TESTS:.nec=.outcpp) $(DO_TESTS:.nec=.outc) $(DO_TESTS:.nec=.diff)
	echo "Test Complete"

############################################################################
#	
#	PYTHON TESTHARNESS
#
#	This tests the python bindings
#
python:
	echo "Put something in here to test the python code"

