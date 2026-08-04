mfpic_end_to_end_test("spatially_constant_maxwellian")

file(
    CREATE_LINK
    ${CMAKE_CURRENT_SOURCE_DIR}/spatially_constant_maxwellian_empire.yaml
    ${CMAKE_CURRENT_BINARY_DIR}/spatially_constant_maxwellian_empire.yaml
    SYMBOLIC)