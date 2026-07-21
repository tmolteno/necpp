#####################################################################################
#
#   CMake build for nec2++ — replaces the hand-written Makefile.
#
#   Copyright (C) 2004-2025  Timothy C.A. Molteno
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#####################################################################################

cmake_minimum_required(VERSION 3.16)

project(necpp
    VERSION 2.1.1
    DESCRIPTION "NEC2++ antenna modelling library and tools"
    HOMEPAGE_URL "https://github.com/tmolteno/necpp"
    LANGUAGES CXX)

# --- C++ standard ----------------------------------------------------------
# C++17 is required by the codebase (structured bindings, etc.).
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to a Release-flavoured build when the caller did not pick one.
# Multi-config generators (Visual Studio, Xcode) choose per-config instead.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Position-independent code so the same objects can back the shared library
# and be linked into executables on platforms that require it.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Build date stamped into the version string written to config.h.
string(TIMESTAMP NECPP_BUILD_DATE "%Y-%m-%d")

# --- Options ---------------------------------------------------------------
option(NECPP_BUILD_TESTS    "Build the Catch2 unit tests"                  ON)
option(NECPP_BUILD_WASM     "Build the Emscripten/WASM target (needs emcc)" OFF)
option(NECPP_BUILD_DOCS     "Build the Doxygen documentation"              OFF)
option(BUILD_SHARED_LIBS    "Also build a shared libnecpp in addition to the static one" ON)

# Bundled Eigen is always used; there is no external-Eigen path.
set(NECPP_EIGEN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/eigen)

# Generate config.h from the template. Output lives in ${CMAKE_BINARY_DIR}
# and is put on the include path by every target that needs it.
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/config.h.in
    ${CMAKE_BINARY_DIR}/config.h
    @ONLY)

add_subdirectory(src)

if(NECPP_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# --- pkg-config file -------------------------------------------------------
# Generated alongside the build so downstream users of `pkg-config --cflags
# --libs necpp` have a path to query. Installed by src/CMakeLists.txt.
include(InstallRequiredSystemLibraries)

# --- CPack -----------------------------------------------------------------
# Replaces the stale debian/ packaging. Produces .deb, .rpm, and .tar.xz.
include(CPackComponent)

set(CPACK_PACKAGE_NAME              "necpp")
set(CPACK_PACKAGE_VENDOR            "Timothy C.A. Molteno")
set(CPACK_PACKAGE_CONTACT           "tim@molteno.net")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "NEC2++ antenna modelling system")
set(CPACK_PACKAGE_DESCRIPTION       "NEC2++ (Numerical Electromagnetics Code) is software for modelling antennas using the Method of Moments. It is a complete rewrite of the original NEC-2 FORTRAN software developed by Lawrence Livermore Laboratories.")
set(CPACK_PACKAGE_HOMEPAGE_URL      "https://github.com/tmolteno/necpp")
set(CPACK_RESOURCE_FILE_LICENSE     ${CMAKE_CURRENT_SOURCE_DIR}/COPYING)
set(CPACK_PACKAGE_VERSION_MAJOR     ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR     ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH     ${PROJECT_VERSION_PATCH})
set(CPACK_PACKAGE_VERSION           ${PROJECT_VERSION})

# Split the install into runtime (lib), development (headers), and CLI
# components so packagers can ship libnecpp / necpp-dev / necpp separately.
# Note: no GROUP here — putting components in a GROUP makes CPack merge them
# into a single group-named package on DEB/RPM, which defeats the split.
cpack_add_component(Runtime
    DISPLAY_NAME "nec2++ runtime library"
    DESCRIPTION  "Shared library libnecpp.so for running programs built against nec2++.")
cpack_add_component(Development
    DISPLAY_NAME "development files"
    DESCRIPTION  "Static library, headers, pkg-config and CMake package files for compiling applications against nec2++.")
cpack_add_component(CLI
    DISPLAY_NAME "nec2++ command-line tools"
    DESCRIPTION  "The nec2++ and nec2diff command-line programs.")

# Debian-specific knobs.
set(CPACK_DEBIAN_PACKAGE_MAINTAINER   ${CPACK_PACKAGE_CONTACT})
set(CPACK_DEBIAN_RUNTIME_PACKAGE_NAME "libnecpp${PROJECT_VERSION_MAJOR}")
set(CPACK_DEBIAN_DEVELOPMENT_PACKAGE_NAME "libnecpp-dev")
set(CPACK_DEBIAN_CLI_PACKAGE_NAME     "necpp")
set(CPACK_DEB_COMPONENT_INSTALL       ON)

# RPM-specific knobs.
set(CPACK_RPM_PACKAGE_LICENSE         "GPL-2.0-or-later")
set(CPACK_RPM_RUNTIME_PACKAGE_NAME    "libnecpp${PROJECT_VERSION_MAJOR}")
set(CPACK_RPM_DEVELOPMENT_PACKAGE_NAME "libnecpp-devel")
set(CPACK_RPM_CLI_PACKAGE_NAME        "necpp")
set(CPACK_RPM_COMPONENT_INSTALL       ON)

include(CPack)
