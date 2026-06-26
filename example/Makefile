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

# Eigen is bundled in ../src/eigen3/ — no external Eigen or LAPACK needed.
# Build with:  make test_cpp

all: test_cpp

CONFIG_H = config.h

$(CONFIG_H):
	echo '#ifndef CONFIG_H' > $@
	echo '#define CONFIG_H' >> $@
	echo '#define VERSION "2.0.0"' >> $@
	echo '#define BUILD_DATE "unknown"' >> $@
	echo '#endif' >> $@

NECPP_SRC  = ../src/c_ggrid.cpp ../src/matrix_algebra.cpp ../src/nec_ground.cpp \
             ../src/c_evlcom.cpp ../src/c_plot_card.cpp ../src/misc.cpp \
             ../src/nec_context.cpp ../src/nec_output.cpp ../src/c_geometry.cpp \
             ../src/nec_exception.cpp ../src/nec_radiation_pattern.cpp \
             ../src/nec_structure_currents.cpp ../src/electromag.cpp ../src/nec_results.cpp \
             ../src/libNEC.cpp
NECPP_OBJS = $(NECPP_SRC:.cpp=.o)

CPP_SRC  = test_cpp.cpp
CPP_OBJS = $(CPP_SRC:.cpp=.o)

CXXFLAGS = -std=c++17 -I ../src/ -I ../src/eigen3 -I . -O2 -Wall

test_cpp: $(CONFIG_H) ${CPP_OBJS} ${NECPP_OBJS}
	$(CXX) -o test_cpp ${CPP_OBJS} ${NECPP_OBJS} -lstdc++ -lm
	./test_cpp

cpp_clean:
	rm -f ${CPP_OBJS} ${NECPP_OBJS}

clean: cpp_clean
	rm -f test_cpp test_c $(CONFIG_H)
