# CMake script run by the necpp_smoke_hertzian_dipole ctest.
# Runs nec2++ -i <input> and fails if "TOTAL RUN TIME" is missing from output.
execute_process(
    COMMAND ${NEC2PP} -i ${INPUT} -o ${CMAKE_CURRENT_BINARY_DIR}/smoke.out
    RESULT_VARIABLE rc
    OUTPUT_VARIABLE stdout
    ERROR_VARIABLE stderr)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "nec2++ exited with code ${rc}\nstderr:\n${stderr}")
endif()
file(READ ${CMAKE_CURRENT_BINARY_DIR}/smoke.out out)
if(NOT out MATCHES "TOTAL RUN TIME")
    message(FATAL_ERROR "nec2++ output did not contain 'TOTAL RUN TIME'")
endif()
message(STATUS "Smoke test passed: ${INPUT}")
