add_executable(mfpic mfpic.cpp)
target_include_directories(mfpic PRIVATE
  ${CMAKE_SOURCE_DIR}
)
target_link_libraries(mfpic
  libmfpic
  mfem
)
install(TARGETS mfpic DESTINATION bin)
