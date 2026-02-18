set(test_prefix "end_to_end_tests:sod_shock_tube")

file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/sod_shock_tube.py ${CMAKE_CURRENT_BINARY_DIR}/sod_shock_tube.py SYMBOLIC)

add_test(
  NAME "${test_prefix}:run"
  COMMAND ${Python3_EXECUTABLE} sod_shock_tube.py run $<TARGET_FILE:mfpic>
)

add_test(
  NAME "${test_prefix}:analyze"
  COMMAND ${Python3_EXECUTABLE} sod_shock_tube.py analyze
)

set_tests_properties("${test_prefix}:run" PROPERTIES FIXTURES_SETUP ${test_prefix})
set_tests_properties("${test_prefix}:analyze" PROPERTIES FIXTURES_REQUIRED ${test_prefix})