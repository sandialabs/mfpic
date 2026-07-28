set(end_to_end_test_num_threads 2)

function(mfpic_end_to_end_test test_name)
  if(mfpic_enable_testing)
    set(test_prefix "end_to_end_tests:${test_name}")

    file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.py ${CMAKE_CURRENT_BINARY_DIR}/${test_name}.py SYMBOLIC)

    add_test(
      NAME "${test_prefix}:run"
      COMMAND ${Python3_EXECUTABLE} ${test_name}.py run $<TARGET_FILE:mfpic>
    )

    add_test(
      NAME "${test_prefix}:analyze"
      COMMAND ${Python3_EXECUTABLE} ${test_name}.py analyze
    )

    set_tests_properties("${test_prefix}:run" PROPERTIES
      FIXTURES_SETUP ${test_prefix}
      PROCESSORS ${end_to_end_test_num_threads}
      ENVIRONMENT "OMP_NUM_THREADS=${end_to_end_test_num_threads}"
    )
    set_tests_properties("${test_prefix}:analyze" PROPERTIES FIXTURES_REQUIRED ${test_prefix})
  endif()
endfunction()

file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/python ${CMAKE_CURRENT_BINARY_DIR}/python SYMBOLIC)

add_subdirectory(advecting_gaussian)
add_subdirectory(charged_shock_tube)
add_subdirectory(expanding_gaussian)
add_subdirectory(sheath)
add_subdirectory(sod_shock_tube)
add_subdirectory(twostream_instability)
add_subdirectory(twostream_instability_fluid)
