file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/read_mesh_data.py ${CMAKE_CURRENT_BINARY_DIR}/read_mesh_data.py SYMBOLIC)

add_subdirectory(sheath)
add_subdirectory(sod_shock_tube)
add_subdirectory(twostream_instability)
