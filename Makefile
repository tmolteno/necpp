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
EIGEN_DIR = $(SRC_DIR)/eigen3
BUILD_DIR = build/simple

VERSION   = 2.1.1
BUILD_DATE = $(shell date +"%Y-%m-%d")

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
	@echo '#define VERSION "$(VERSION)"'       >> $@
	@echo '#define BUILD_DATE "$(BUILD_DATE)"' >> $@
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

# --- Default targets ---

clean:
	rm -rf $(BUILD_DIR) nec2++ nec2diff

install: nec2++
	install -d $(DESTDIR)/usr/local/bin
	install nec2++ $(DESTDIR)/usr/local/bin/

.PHONY: all clean install
