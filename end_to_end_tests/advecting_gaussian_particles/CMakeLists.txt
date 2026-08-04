mfpic_end_to_end_test("advecting_gaussian_particles")

file(
    CREATE_LINK
    ${CMAKE_CURRENT_SOURCE_DIR}/compute_velocity_histogram.py
    ${CMAKE_CURRENT_BINARY_DIR}/compute_velocity_histogram.py
    SYMBOLIC)