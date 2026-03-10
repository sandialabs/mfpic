set(test_prefix "end_to_end_tests:twostream_instability_fluid")

file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/twostream_instability_fluid.py ${CMAKE_CURRENT_BINARY_DIR}/twostream_instability_fluid.py SYMBOLIC)

add_test(
  NAME "${test_prefix}:run"
  COMMAND ${Python3_EXECUTABLE} twostream_instability_fluid.py run $<TARGET_FILE:mfpic>
)

add_test(
  NAME "${test_prefix}:analyze"
  COMMAND ${Python3_EXECUTABLE} twostream_instability_fluid.py analyze
)

set_tests_properties("${test_prefix}:run" PROPERTIES FIXTURES_SETUP ${test_prefix})
set_tests_properties("${test_prefix}:analyze" PROPERTIES FIXTURES_REQUIRED ${test_prefix})
