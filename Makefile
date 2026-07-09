#####################################################################################
#
#   Simple Makefile for nec2++ — no autotools required
#
#   Copyright (C) 2004-2025  Timothy C.A. Molteno
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#####################################################################################

# Default target (must be first)
all: nec2++ nec2diff

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wshadow
LDFLAGS  ?= -lm

# Add bounds checking:  make DEBUG=1
ifdef DEBUG
  CXXFLAGS += -O0 -g3 -DNEC_ERROR_CHECK=1
endif

# Add typesafe integer checking:  make TYPECHECK=1
ifdef TYPECHECK
  CXXFLAGS += -DTYPESAFE_PEDANTIC=1
endif

SRC_DIR   = src
EIGEN_DIR = $(SRC_DIR)/eigen
BUILD_DIR = build/simple

NECPP_VERSION   = 2.1.1
NECPP_BUILD_DATE = $(shell date +"%Y-%m-%d")

INCLUDES  = -I $(SRC_DIR) -isystem $(EIGEN_DIR) -I $(BUILD_DIR)

LIB_SRCS  = c_evlcom c_geometry c_ggrid c_plot_card electromag libNEC \
            matrix_algebra misc nec_context nec_exception nec_ground \
            nec_output nec_radiation_pattern nec_results nec_structure_currents

LIB_OBJS  = $(addprefix $(BUILD_DIR)/, $(addsuffix .o, $(LIB_SRCS)))

# --- Build directory ---
$(BUILD_DIR):
	@mkdir -p $@

# --- config.h ---
$(BUILD_DIR)/config.h: | $(BUILD_DIR)
	@echo '#ifndef CONFIG_H'           >  $@
	@echo '#define CONFIG_H'           >> $@
	@echo '#define NECPP_VERSION "$(NECPP_VERSION)"'       >> $@
	@echo '#define NECPP_BUILD_DATE "$(NECPP_BUILD_DATE)"' >> $@
	@echo '#endif'                     >> $@

# --- Compile .cpp -> .o ---
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(BUILD_DIR)/config.h | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# --- nec2++ binary ---
nec2++: $(LIB_OBJS) $(SRC_DIR)/nec2cpp.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ \
	  $(SRC_DIR)/nec2cpp.cpp $(SRC_DIR)/XGetopt.cpp $(LIB_OBJS) $(LDFLAGS) -lstdc++

# --- nec2diff utility ---
nec2diff: $(SRC_DIR)/necDiff.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $(SRC_DIR)/necDiff.cpp $(LDFLAGS) -lstdc++

# --- Test runner ---
TEST_OBJS  = safe_array_tb matrix_algebra_tb c_geometry_tb c_evlcom_tb nec_context_tb nec2cpp_tb
TEST_CXXFLAGS = -std=c++17 -I $(SRC_DIR) -isystem $(EIGEN_DIR) -I $(BUILD_DIR) -DNEC_ERROR_CHECK=1 -g

test: $(LIB_OBJS)
	# Recompile library with bounds checking
	@mkdir -p $(BUILD_DIR)
	@for src in $(LIB_SRCS) XGetopt; do \
		$(CXX) $(TEST_CXXFLAGS) -c $(SRC_DIR)/$$src.cpp -o $(BUILD_DIR)/test_$$src.o; \
	done
	# Compile nec2cpp with renamed main
	$(CXX) $(TEST_CXXFLAGS) -Dmain=nec2cpp_main_renamed -c $(SRC_DIR)/nec2cpp.cpp -o $(BUILD_DIR)/test_nec2cpp.o
	# Compile test files
	@for tb in $(TEST_OBJS); do \
		$(CXX) $(TEST_CXXFLAGS) -c $(SRC_DIR)/$$tb.cpp -o $(BUILD_DIR)/$$tb.o; \
	done
	# Build Catch2 main
	@printf '#define CATCH_CONFIG_MAIN\n#include "catch.hpp"\n' > $(BUILD_DIR)/test_main.cpp
	$(CXX) $(TEST_CXXFLAGS) -c $(BUILD_DIR)/test_main.cpp -o $(BUILD_DIR)/test_main.o
	# Link and run
	$(CXX) -o $(BUILD_DIR)/test_runner $(BUILD_DIR)/test_main.o \
		$(addprefix $(BUILD_DIR)/, $(TEST_OBJS:=.o)) \
		$(BUILD_DIR)/test_nec2cpp.o \
		$(addprefix $(BUILD_DIR)/test_, $(LIB_SRCS:=.o)) \
		$(BUILD_DIR)/test_XGetopt.o \
		-Wl,--allow-multiple-definition -lm -lstdc++
	./$(BUILD_DIR)/test_runner "[safe_array],[lu_decompose],[factrs_solves],[bessel],[card_handler],[compat],[readmn],[load_line],[parser_dispatch],[helix],[intersection_check],[false_intersection],[parse_nec_card]" -s

# --- Default targets ---

clean:
	rm -rf $(BUILD_DIR) nec2++ nec2diff

install: nec2++
	install -d $(DESTDIR)/usr/local/bin
	install nec2++ $(DESTDIR)/usr/local/bin/

# --- WASM target (requires Emscripten: emcc) ---
# Build with:  make wasm EMCC=emcc
EMCC ?= emcc
WASM_EXPORTS = '["_nec_create_context","_nec_delete_context","_nec_process_input","_nec_get_output","_nec_get_output_length","_nec_free"]'
WASM_RUNTIME  = '["ccall","cwrap","UTF8ToString","lengthBytesUTF8"]'

wasm: $(LIB_OBJS) $(SRC_DIR)/nec_wasm.cpp
	$(EMCC) -std=c++17 -O2 $(INCLUDES) \
		-s WASM=1 \
		-s EXPORTED_FUNCTIONS=$(WASM_EXPORTS) \
		-s EXPORTED_RUNTIME_METHODS=$(WASM_RUNTIME) \
		-s ALLOW_MEMORY_GROWTH=1 \
		--no-entry \
		-o nec2pp.js \
		$(SRC_DIR)/nec_wasm.cpp $(SRC_DIR)/XGetopt.cpp $(LIB_OBJS) -lstdc++
	@echo "WASM build complete: nec2pp.js + nec2pp.wasm"

.PHONY: all clean install test wasm
