#####################################################################################
#
#	Makefile showing how to use libnecpp, the library
#	version of nec++. This lets you incorporate NEC2
#	antenna electromagnetic modeling into your own programs
#
#	Copyright (C) 2004-2005, 2015  Timothy C.A. Molteno
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

LD_LIBRARY_PATH=/usr/lib/atlas-base/atlas solve
# TEST_LIBRARY_PATH=../lib
# ./configure  --prefix=`pwd` --with-bounds

all: test_cpp

test_static: test_nec.c
	gcc --static -o test_static test_nec.c -lnecpp -lstdc++ -lm -llapack  -lblas

test_c: test_nec.o
	gcc -o test_c -Wl,-rpath -Wl,$(TEST_LIBRARY_PATH) test_nec.o -lnecpp

test: test_c
	./test_c
	ldd test_c


#####################################################################################
#
#	The target below shows how to compile and use the nec2++ codebase
#	from within other C++ programs (GPL compatable of course!)
#

NECPP_SRC  = ../src/c_ggrid.cpp ../src/matrix_algebra.cpp ../src/nec_ground.cpp \
             ../src/c_evlcom.cpp ../src/c_plot_card.cpp ../src/misc.cpp \
             ../src/nec_context.cpp ../src/nec_output.cpp ../src/c_geometry.cpp \
             ../src/nec_exception.cpp ../src/nec_radiation_pattern.cpp \
	     ../src/nec_structure_currents.cpp ../src/electromag.cpp ../src/nec_results.cpp
NECPP_OBJS = $(NECPP_SRC:.cpp=.o)


CPP_SRC  = test_cpp.cpp
CPP_OBJS = $(CPP_SRC:.cpp=.o)

CXXFLAGS=-I ../src/ -I .. -O2 -Wall


test_cpp: ${CPP_OBJS} ${NECPP_OBJS}
	gcc -o test_cpp ${CPP_OBJS} ${NECPP_OBJS} -lstdc++ -lm -llapack
	./test_cpp

cpp_clean:
	rm -f ${CPP_OBJS} ${NECPP_OBJS} 

clean: cpp_clean
	rm -f test_cpp test_static test_c
