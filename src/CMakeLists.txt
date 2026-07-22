#####################################################################################
#
#   Library, executables, install/export for nec2++.
#
#####################################################################################

# --- The canonical nec2++ source list. -------------------------------------
# Single source of truth: every target below (static lib, shared lib, tests,
# WASM) is built from this one list. Add/remove a .cpp here and it propagates
# everywhere.
set(NECPP_LIB_SRCS
    c_evlcom.cpp
    c_geometry.cpp
    c_ggrid.cpp
    c_plot_card.cpp
    electromag.cpp
    libNEC.cpp
    matrix_algebra.cpp
    misc.cpp
    nec_context.cpp
    nec_exception.cpp
    nec_ground.cpp
    nec_output.cpp
    nec_radiation_pattern.cpp
    nec_results.cpp
    nec_structure_currents.cpp
)

# Public-facing headers installed alongside the library.
set(NECPP_PUBLIC_HEADERS
    libnecpp.h
    nec_context.h
    c_geometry.h
    nec_ground.h
    nec_radiation_pattern.h
    nec_results.h
    nec_structure_currents.h
    nec_output.h
    nec_exception.h
    c_evlcom.h
    c_ggrid.h
    c_plot_card.h
    electromag.h
    matrix_algebra.h
    misc.h
    nec_wire.h
    math_util.h
    safe_array.h
    common.h
    typesafe_stdint.h
    nec_debug.h
    nec_card_parser.h
    AntennaInput.h
    BaseInput.h
    CurrentInput.h
    PowerBudget.h
    RadiationInput.h
    nec2cpp.h
    XGetopt.h
)

# --- OBJECT library: compiled once, linked into both static and shared. ----
add_library(necpp_obj OBJECT ${NECPP_LIB_SRCS})
target_include_directories(necpp_obj
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${NECPP_EIGEN_DIR}
        ${CMAKE_BINARY_DIR})        # for generated config.h
# Warning flags. The old Makefile used -Wall -Wextra -Wshadow (GCC/Clang).
# MSVC uses /W3 or /W4 instead; enable /W3 as a reasonable baseline.
if(MSVC)
    target_compile_options(necpp_obj PRIVATE /W3)
else()
    target_compile_options(necpp_obj PRIVATE -Wall -Wextra -Wshadow)
endif()

# --- Static library --------------------------------------------------------
add_library(necpp_static STATIC $<TARGET_OBJECTS:necpp_obj>)
set_target_properties(necpp_static PROPERTIES
    OUTPUT_NAME necpp
    EXPORT_NAME necpp_static
    POSITION_INDEPENDENT_CODE ON)
target_include_directories(necpp_static
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${NECPP_EIGEN_DIR}>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
        $<INSTALL_INTERFACE:include/necpp>)
# On Linux, link the math library.
target_link_libraries(necpp_static PUBLIC ${CMAKE_DL_LIBS})
if(UNIX AND NOT APPLE)
    target_link_libraries(necpp_static PUBLIC m)
endif()

# --- Shared library (optional) ---------------------------------------------
if(BUILD_SHARED_LIBS)
    add_library(necpp_shared SHARED $<TARGET_OBJECTS:necpp_obj>)
    set_target_properties(necpp_shared PROPERTIES
        OUTPUT_NAME necpp
        EXPORT_NAME necpp
        SOVERSION ${PROJECT_VERSION_MAJOR}
        VERSION ${PROJECT_VERSION})
    target_include_directories(necpp_shared
        PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
            $<BUILD_INTERFACE:${NECPP_EIGEN_DIR}>
            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
            $<INSTALL_INTERFACE:include/necpp>)
    target_link_libraries(necpp_shared PUBLIC ${CMAKE_DL_LIBS})
    if(UNIX AND NOT APPLE)
        target_link_libraries(necpp_shared PUBLIC m)
    endif()
    add_library(necpp::necpp ALIAS necpp_shared)
    # Default library target consumed by executables: prefer the shared lib
    # when it exists, otherwise the static one.
    set(NECPP_LIB_TARGET necpp_shared)
else()
    # Static-only build: export the static target as the canonical necpp::necpp
    # so downstream find_package(necpp) users get one stable name either way.
    set_target_properties(necpp_static PROPERTIES EXPORT_NAME necpp)
    add_library(necpp::necpp ALIAS necpp_static)
    set(NECPP_LIB_TARGET necpp_static)
endif()

# --- Command-line executables ----------------------------------------------
# nec2++ — the main NEC-2 simulation driver.
add_executable(nec2++ nec2cpp.cpp XGetopt.cpp)
target_link_libraries(nec2++ PRIVATE ${NECPP_LIB_TARGET})
set_target_properties(nec2++ PROPERTIES OUTPUT_NAME "nec2++")

# nec2diff — utility for diffing two NEC output files.
add_executable(nec2diff necDiff.cpp)
target_link_libraries(nec2diff PRIVATE ${NECPP_LIB_TARGET})

# Under a WASM build we only want the Emscripten artifacts; skip the native
# executables so `cmake --build` doesn't waste time compiling nec2++/nec2diff
# as (broken) .js shells.
if(NECPP_BUILD_WASM)
    set_target_properties(nec2++ nec2diff PROPERTIES EXCLUDE_FROM_ALL ON)
endif()

# --- WASM target (opt-in via -DNECPP_BUILD_WASM=ON, needs emcc) ------------
if(NECPP_BUILD_WASM)
    add_executable(nec2pp_wasm
        ${CMAKE_CURRENT_SOURCE_DIR}/nec_wasm.cpp
        XGetopt.cpp
        ${NECPP_LIB_SRCS})
    target_include_directories(nec2pp_wasm PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${NECPP_EIGEN_DIR}
        ${CMAKE_BINARY_DIR})
    # Flags carried verbatim from the old Makefile wasm target.
    target_link_options(nec2pp_wasm PRIVATE
        -sWASM=1
        "-sEXPORTED_FUNCTIONS=[\"_nec_create_context\",\"_nec_delete_context\",\"_nec_process_input\",\"_nec_get_output\",\"_nec_get_output_length\",\"_nec_free\"]"
        "-sEXPORTED_RUNTIME_METHODS=[\"ccall\",\"cwrap\",\"UTF8ToString\",\"lengthBytesUTF8\"]"
        -sALLOW_MEMORY_GROWTH=1
        --no-entry)
    # OUTPUT_NAME nec2pp + SUFFIX .js → Emscripten emits nec2pp.js + nec2pp.wasm,
    # matching the legacy Makefile output names.
    set_target_properties(nec2pp_wasm PROPERTIES
        OUTPUT_NAME nec2pp
        PREFIX ""
        SUFFIX ".js")
    # Surface the wasm next to the js (Emscripten emits it in the build dir).
    add_custom_command(TARGET nec2pp_wasm POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_CURRENT_BINARY_DIR}/nec2pp.wasm
            ${CMAKE_BINARY_DIR}/nec2pp.wasm
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_CURRENT_BINARY_DIR}/nec2pp.js
            ${CMAKE_BINARY_DIR}/nec2pp.js
        COMMENT "Copying WASM artifacts to build root")
endif()

# --- Installation ----------------------------------------------------------
include(GNUInstallDirs)

install(TARGETS necpp_static
    EXPORT necppTargets
    COMPONENT Development
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

if(BUILD_SHARED_LIBS)
    install(TARGETS necpp_shared
        EXPORT necppTargets
        COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

install(TARGETS nec2++ nec2diff
    COMPONENT CLI
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Public headers → include/necpp/
install(FILES
    ${NECPP_PUBLIC_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/necpp
    COMPONENT Development)

# Bundled Eigen headers needed to use the library (it is header-only and
# types from it leak into the public API).
install(DIRECTORY ${NECPP_EIGEN_DIR}/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/necpp/eigen
    COMPONENT Development
    FILES_MATCHING PATTERN "*.h"
                    PATTERN "*.hpp"
                    PATTERN "CMakeLists.txt" EXCLUDE)

# Man pages (exist in the repo root, were never installed by the old Makefile).
install(FILES
    ${CMAKE_SOURCE_DIR}/nec2++.1
    ${CMAKE_SOURCE_DIR}/nec2diff.1
    DESTINATION ${CMAKE_INSTALL_MANDIR}/man1
    COMPONENT CLI OPTIONAL)

# --- pkg-config file -------------------------------------------------------
configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/necpp.pc.in
    ${CMAKE_BINARY_DIR}/necpp.pc
    @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/necpp.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
    COMPONENT Development)

# --- CMake package config (find_package(necpp)) ----------------------------
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_BINARY_DIR}/necppConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)
configure_package_config_file(
    ${CMAKE_SOURCE_DIR}/cmake/necppConfig.cmake.in
    ${CMAKE_BINARY_DIR}/necppConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/necpp)

install(EXPORT necppTargets
    FILE necppTargets.cmake
    NAMESPACE necpp::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/necpp
    COMPONENT Development)

install(FILES
    ${CMAKE_BINARY_DIR}/necppConfig.cmake
    ${CMAKE_BINARY_DIR}/necppConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/necpp
    COMPONENT Development)
