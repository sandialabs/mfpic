cmake_minimum_required(VERSION 3.19)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(mfpic)

if(MSVC)
  add_compile_options(/W4 /WX)
else()
  add_compile_options(-Wall -Wextra -Wpedantic -Werror)
endif()

find_package(HDF5)

find_package(GTest)
find_package(Python3 COMPONENTS Interpreter)
set(mfpic_enable_testing False)
if(GTest_FOUND AND Python3_FOUND)
  set(mfpic_enable_testing True)
  enable_testing()
  include(GoogleTest)
  add_subdirectory(end_to_end_tests)
endif()
function(mfpic_unit_test test_name)
  if(mfpic_enable_testing)
    add_executable("${test_name}.test" "${test_name}.test.cpp")
    target_include_directories("${test_name}.test" PRIVATE
      ${CMAKE_SOURCE_DIR}
    )
    target_include_directories("${test_name}.test" SYSTEM PRIVATE
      mfem
    )
    target_link_libraries("${test_name}.test"
      libmfpic
      mfem
      GTest::gtest_main
    )
    gtest_discover_tests("${test_name}.test")
  endif()
endfunction()

add_subdirectory(libmfpic)
add_subdirectory(mfpic)
