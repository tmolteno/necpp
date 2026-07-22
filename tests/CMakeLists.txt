#####################################################################################
#
#   Catch2 v3 unit tests for nec2++.
#
#   Catch2 is obtained via FetchContent at configure time (no vendored header).
#   The library sources are compiled with NEC_ERROR_CHECK=1, and nec2cpp.cpp
#   is linked with its main() renamed so Catch owns main().
#
#####################################################################################

include(FetchContent)

# Catch2 v3 — downloaded at configure time, provides the test framework
# and main() via Catch2::Catch2WithMain.
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.7.1
    GIT_SHALLOW    TRUE
)
# Suppress the "Catch2" project's find_package interference.
set(CATCH_INSTALL_DOCS OFF CACHE BOOL "" FORCE)
set(CATCH_INSTALL_HELPERS OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(Catch2)

# All 7 test files (including math_util_tb.cpp, which the old Makefile
# excluded for historical reasons — it has a legitimate [nec_3vector] test).
set(NECPP_TEST_SRCS
    ${CMAKE_SOURCE_DIR}/src/safe_array_tb.cpp
    ${CMAKE_SOURCE_DIR}/src/matrix_algebra_tb.cpp
    ${CMAKE_SOURCE_DIR}/src/c_geometry_tb.cpp
    ${CMAKE_SOURCE_DIR}/src/c_evlcom_tb.cpp
    ${CMAKE_SOURCE_DIR}/src/math_util_tb.cpp
    ${CMAKE_SOURCE_DIR}/src/nec_context_tb.cpp
    ${CMAKE_SOURCE_DIR}/src/nec2cpp_tb.cpp
)

# Test runner = Catch2 main + test sources + nec2cpp.cpp (main renamed) + the
# library sources, all compiled under the bounds-checking defines. This matches
# the old Makefile, which recompiled the whole library with -DNEC_ERROR_CHECK=1
# rather than linking the production objects.
add_executable(nec2++_tests
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
    ${CMAKE_BINARY_DIR})

# Bounds checking on (matches the old Makefile TEST_CXXFLAGS).
target_compile_definitions(nec2++_tests PRIVATE NEC_ERROR_CHECK=1)

# Rename nec2cpp's main() so the Catch runner owns main().
set_source_files_properties(${CMAKE_SOURCE_DIR}/src/nec2cpp.cpp PROPERTIES
    COMPILE_DEFINITIONS "main=nec2cpp_main_renamed")

target_link_libraries(nec2++_tests PRIVATE
    Catch2::Catch2WithMain)
if(UNIX AND NOT APPLE)
    target_link_libraries(nec2++_tests PRIVATE m)
endif()

# Run all test cases except [surface_patch], which is a pre-existing failure
# (excluded from the old CTest tag filter since before the CMake migration).
# The ~tag syntax is Catch2 v3's tag exclusion — all future test cases are
# picked up automatically without maintaining a whitelist.
add_test(NAME necpp_unit
    COMMAND nec2++_tests "~[surface_patch]" -s)

# Smoke test: run a real simulation through the binary and check it produced
# the expected output marker.
add_test(NAME necpp_smoke_hertzian_dipole
    COMMAND ${CMAKE_COMMAND}
        -DNEC2PP=$<TARGET_FILE:nec2++>
        -DINPUT=${CMAKE_SOURCE_DIR}/testharness/data/herzian_dipole.nec
        -P ${CMAKE_CURRENT_SOURCE_DIR}/smoke_test.cmake)