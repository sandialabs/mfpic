file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/python ${CMAKE_CURRENT_BINARY_DIR}/python SYMBOLIC)

add_subdirectory(sheath)
add_subdirectory(sod_shock_tube)
add_subdirectory(twostream_instability)
add_subdirectory(twostream_instability_fluid)
