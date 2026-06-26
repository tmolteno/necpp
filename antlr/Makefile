########################################
#
#	Build the NEC2++ parser using the
#	ANTLR compiler construction toolkit
#
#	DEPRECATED: This is experimental code.
#	The production parser is in src/nec2cpp.cpp.
#	This Makefile is kept for reference only.
#

.SUFFIXES: .parse .nec

NECPP_SRC 	= ../src/c_evlcom.cpp ../src/c_geometry.cpp ../src/c_ggrid.cpp \
		../src/c_plot_card.cpp ../src/matrix_algebra.cpp ../src/misc.cpp \
		../src/nec_context.cpp ../src/nec_ground.cpp ../src/nec_output.cpp \
		../src/nec_radiation_pattern.cpp ../src/nec_exception.cpp \
		../src/electromag.cpp ../src/nec_structure_currents.cpp ../src/nec_results.cpp \
		../src/XGetopt.cpp ../src/libNEC.cpp
NECPP_OBJ 	= $(NECPP_SRC:.cpp=.o)

ANTLR_SRC 	= antlr_test.cpp NECLexer.cpp NECParser.cpp
ANTLR_OBJ 	= $(ANTLR_SRC:.cpp=.o)
ANTLR		= java -classpath /usr/share/java/antlr.jar antlr.Tool

all:	antlr test

echo:
	echo ${NECPP_SRC}

antlr:	nec.g
	${ANTLR} nec.g

clean:
	rm -f ${ANTLR_OBJ} ${NECPP_OBJ} NEC* antlr_test

LDFLAGS = -lstdc++ -lantlr -lm
CXXFLAGS = -std=c++17 -I ../ -I ../src/ -I ../src/eigen3 -g3 -O0 -fPIC

test: ${ANTLR_OBJ} ${ANTLR_SRC} ${NECPP_OBJ} ${NECPP_SRC}
	g++ -o antlr_test ${ANTLR_OBJ} ${NECPP_OBJ} ${LDFLAGS}
	- ./antlr_test -i ../testharness/data/example3.nec -o example3.out

TEST_FILES := $(shell find ../testharness/data/ -name "*.nec")
TEST_OUT	= $(TEST_FILES:.nec=.antlr)

%.nec %.antlr:
	- ./antlr_test -i $*.nec -o $@

test_clean:
	rm -f ${TEST_OUT}

test_all: ${TEST_OUT}
	echo "DONE"
