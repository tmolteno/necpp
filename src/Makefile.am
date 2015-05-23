# 
#   Copyright (C) 2004,2015  Timothy C.A. Molteno
#   
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#   
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
bin_PROGRAMS = nec2diff nec2++

nec2___SOURCES = common.h nec2cpp.cpp nec2cpp.h XGetopt.cpp XGetopt.h

nec2diff_SOURCES =  AntennaInput.h BaseInput.h CurrentInput.h necDiff.cpp PowerBudget.h RadiationInput.h


include_HEADERS = libnecpp.h 
noinst_HEADERS = nec_debug.h nec_exception.h nec_wire.h

# set the include path found by configure
#INCLUDES = $(all_includes)
AM_CPPFLAGS = $(all_includes)

# the library search path.
nec2___LDADD = libnecpp.la

lib_LTLIBRARIES = libnecpp.la

libnecpp_la_SOURCES = c_evlcom.cpp c_geometry.cpp c_ggrid.cpp c_plot_card.cpp \
	libNEC.cpp matrix_algebra.cpp misc.cpp nec_context.cpp nec_exception.cpp \
	nec_ground.cpp nec_output.cpp nec_radiation_pattern.cpp c_evlcom.h c_geometry.h \
	c_ggrid.h c_plot_card.h common.h electromag.h math_util.h matrix_algebra.h misc.h \
	libnecpp.h nec_context.h nec_exception.h nec_ground.h nec_output.h \
	nec_radiation_pattern.h nec_results.h safe_array.h nec_structure_currents.cpp \
	nec_structure_currents.h electromag.cpp nec_results.cpp 
libnecpp_la_CXXFLAGS = -O2 -g0 -Wall -Wextra -Wshadow


libnecpp_la_LDFLAGS = @LIBLAPACK@ -lm -lstdc++
nec2___LDFLAGS = -all-static

# Unit tests
noinst_PROGRAMS = test_manager

test_manager_SOURCES  = $(libnecpp_la_SOURCES) safe_array_tb.cpp matrix_algebra_tb.cpp c_geometry_tb.cpp c_evlcom_tb.cpp
test_manager_CXXFLAGS = -std=c++11 -O0 -g3 -Wall -Wextra -Werror -DNEC_ERROR_CHECK=1
test_manager_LDADD    = $(libnecpp_la_LIBADD) 
test_manager_LDFLAGS    = $(libnecpp_la_LDFLAGS) 




