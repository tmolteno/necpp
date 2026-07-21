#####################################################################################
#
#   Catch2 unit tests for nec2++.
#
#   Mirrors the old Makefile `test` target: the library sources are recompiled
#   with NEC_ERROR_CHECK=1, nec2cpp.cpp is pulled in with its main() renamed
#   so the Catch runner owns main(), and the test sources are linked together.
#
#####################################################################################

# The six test files the old Makefile compiled (see Makefile TEST_OBJS).
# math_util_tb.cpp exists in src/ but was never part of the test runner —
# leave it out to keep behaviour identical.
set(NECPP_TEST_SRCS
    ${CMAKE_SOURCE_DIR}/src/safe_array_tb.cpp
    ${CMAKE_SOURCE_DIR}/src/matrix_algebra_tb.cpp
    ${CMAKE_SOURCE_DIR}/src/c_geometry_tb.cpp
    ${CMAKE_SOURCE_DIR}/src/c_evlcom_tb.cpp
    ${CMAKE_SOURCE_DIR}/src/nec_context_tb.cpp
    ${CMAKE_SOURCE_DIR}/src/nec2cpp_tb.cpp
)

# Test runner = Catch main + test sources + nec2cpp.cpp (main renamed) + the
# library sources, all compiled under the bounds-checking defines. This matches
# the old Makefile, which recompiled the whole library with -DNEC_ERROR_CHECK=1
# rather than linking the production objects.
add_executable(nec2++_tests
    test_main.cpp
    ${NECPP_TEST_SRCS}
    ${CMAKE_SOURCE_DIR}/src/nec2cpp.cpp
    ${CMAKE_SOURCE_DIR}/src/XGetopt.cpp
    ${CMAKE_SOURCE_DIR}/src/c_evlcom.cpp
    ${CMAKE_SOURCE_DIR}/src/c_geometry.cpp
    ${CMAKE_SOURCE_DIR}/src/c_ggrid.cpp
    ${CMAKE_SOURCE_DIR}/src/c_plot_card.cpp
    ${CMAKE_SOURCE_DIR}/src/electromag.cpp
    ${CMAKE_SOURCE_DIR}/src/libNEC.cpp
    ${CMAKE_SOURCE_DIR}/src/matrix_algebra.cpp
    ${CMAKE_SOURCE_DIR}/src/misc.cpp
    ${CMAKE_SOURCE_DIR}/src/nec_context.cpp
    ${CMAKE_SOURCE_DIR}/src/nec_exception.cpp
    ${CMAKE_SOURCE_DIR}/src/nec_ground.cpp
    ${CMAKE_SOURCE_DIR}/src/nec_output.cpp
    ${CMAKE_SOURCE_DIR}/src/nec_radiation_pattern.cpp
    ${CMAKE_SOURCE_DIR}/src/nec_results.cpp
    ${CMAKE_SOURCE_DIR}/src/nec_structure_currents.cpp
)

target_include_directories(nec2++_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${NECPP_EIGEN_DIR}
    ${CMAKE_BINARY_DIR})             # generated config.h + catch.hpp is in src/

# Bounds checking on (matches the old Makefile TEST_CXXFLAGS).
target_compile_definitions(nec2++_tests PRIVATE NEC_ERROR_CHECK=1)

# Rename nec2cpp's main() so the Catch runner (in test_main.cpp) owns main().
# This define must apply ONLY to nec2cpp.cpp — if it leaks into test_main.cpp
# it renames Catch's own main() and the link fails with "undefined reference
# to main". The old Makefile got this right by compiling test_main.cpp in a
# separate step without the define; here we scope it with a source property.
set_source_files_properties(${CMAKE_SOURCE_DIR}/src/nec2cpp.cpp PROPERTIES
    COMPILE_DEFINITIONS "main=nec2cpp_main_renamed")

# nec2cpp.cpp defines objects that also exist in the explicit library sources
# when main is renamed — allow the duplicates just like the old Makefile did.
target_link_options(nec2++_tests PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wl,--allow-multiple-definition>)

target_link_libraries(nec2++_tests PRIVATE m)

# The same Catch tag filter the old Makefile passed on the command line.
set(NECPP_TEST_TAGS
    "[safe_array],[lu_decompose],[factrs_solves],[bessel],[card_handler],\
[compat],[readmn],[load_line],[parser_dispatch],[helix],\
[intersection_check],[false_intersection],[parse_nec_card]")

add_test(NAME necpp_unit
    COMMAND nec2++_tests ${NECPP_TEST_TAGS} -s)

# Smoke test: run a real simulation through the binary and check it produced
# the expected output marker.
add_test(NAME necpp_smoke_hertzian_dipole
    COMMAND ${CMAKE_COMMAND}
        -DNEC2PP=$<TARGET_FILE:nec2++>
        -DINPUT=${CMAKE_SOURCE_DIR}/testharness/data/herzian_dipole.nec
        -P ${CMAKE_CURRENT_SOURCE_DIR}/smoke_test.cmake)
