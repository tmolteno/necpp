bin_PROGRAMS = nec2diff nec2++

nec2___SOURCES = common.h nec2cpp.cpp nec2cpp.h XGetopt.cpp XGetopt.h

nec2diff_SOURCES =  AntennaInput.h BaseInput.h CurrentInput.h necDiff.cpp PowerBudget.h RadiationInput.h

include_HEADERS = libnecpp.h 
noinst_HEADERS = nec_debug.h nec_exception.h nec_wire.h

# set the include path found by configure
INCLUDES = $(all_includes)

# the library search path.
nec2___LDADD = $(top_builddir)/src/libnecpp.la

lib_LTLIBRARIES = libnecpp.la

libnecpp_la_SOURCES = c_evlcom.cpp c_geometry.cpp c_ggrid.cpp c_plot_card.cpp \
	libNEC.cpp matrix_algebra.cpp misc.cpp nec_context.cpp nec_exception.cpp \
	nec_ground.cpp nec_output.cpp nec_radiation_pattern.cpp c_evlcom.h c_geometry.h \
	c_ggrid.h c_plot_card.h common.h electromag.h math_util.h matrix_algebra.h misc.h \
	libnecpp.h nec_context.h nec_exception.h nec_ground.h nec_output.h \
	nec_radiation_pattern.h nec_results.h safe_array.h nec_structure_currents.cpp \
	nec_structure_currents.h electromag.cpp

libnecpp_la_LDFLAGS = @LIBLAPACK@ -lm -lstdc++
nec2___LDFLAGS = -all-static

