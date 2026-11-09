# Regression test for github issue #138: nec2diff treated any two negative
# values as equal (the zero-test compared signed values), reported `inf` for
# pure sign flips, always exited 0, and aborted when called with fewer than
# two file arguments.
#
# The fixtures in data/issue138 are one NEC output file and four single-value
# edits of it: a negative reactance change (-179.39 -> -1.00 ohm), a negative
# dB gain change (-1.71 -> -9.71 dBi), a positive control (49.889 -> 59.889
# ohm) and a reactance sign flip (-179.39 -> +179.39 ohm).
#
# Run as:
#   cmake -DNEC2DIFF=<path> -DDATA_DIR=<dir> -P nec2diff_test.cmake

set(BASE ${DATA_DIR}/issue138_base.out)

# Every edited fixture must be reported as different: exit status 1.
set(DIFFER_CASES
    issue138_reactance
    issue138_gain
    issue138_resistance
    issue138_signflip)

foreach(case_name IN LISTS DIFFER_CASES)
    execute_process(
        COMMAND ${NEC2DIFF} ${BASE} ${DATA_DIR}/${case_name}.out
        OUTPUT_VARIABLE nec2diff_stdout
        ERROR_VARIABLE nec2diff_stderr
        RESULT_VARIABLE nec2diff_result)
    if(NOT nec2diff_result EQUAL 1)
        message(FATAL_ERROR
            "nec2diff ${case_name}: expected exit status 1 (different), got ${nec2diff_result}\n"
            "stdout:\n${nec2diff_stdout}\nstderr:\n${nec2diff_stderr}")
    endif()
endforeach()

# Identical files must compare equal: exit status 0.
execute_process(
    COMMAND ${NEC2DIFF} ${BASE} ${BASE}
    OUTPUT_QUIET ERROR_QUIET
    RESULT_VARIABLE nec2diff_result)
if(NOT nec2diff_result EQUAL 0)
    message(FATAL_ERROR "nec2diff: identical files must exit 0, got ${nec2diff_result}")
endif()

# Fewer than two arguments must print usage and exit 1 (not abort).
execute_process(
    COMMAND ${NEC2DIFF}
    OUTPUT_QUIET ERROR_QUIET
    RESULT_VARIABLE nec2diff_result)
if(NOT nec2diff_result EQUAL 1)
    message(FATAL_ERROR "nec2diff: no arguments must exit 1, got ${nec2diff_result}")
endif()

execute_process(
    COMMAND ${NEC2DIFF} ${BASE}
    OUTPUT_QUIET ERROR_QUIET
    RESULT_VARIABLE nec2diff_result)
if(NOT nec2diff_result EQUAL 1)
    message(FATAL_ERROR "nec2diff: one argument must exit 1, got ${nec2diff_result}")
endif()

message(STATUS "nec2diff issue #138 regression test passed")
