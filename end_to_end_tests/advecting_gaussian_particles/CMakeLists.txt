mfpic_end_to_end_test("advecting_gaussian_particles")

file(
    CREATE_LINK
    ${CMAKE_CURRENT_SOURCE_DIR}/compute_particle_distribution_function.py
    ${CMAKE_CURRENT_BINARY_DIR}/compute_particle_distribution_function.py
    SYMBOLIC)

file(
    CREATE_LINK
    ${CMAKE_CURRENT_SOURCE_DIR}/plot_particle_distribution_function.py
    ${CMAKE_CURRENT_BINARY_DIR}/plot_particle_distribution_function.py
    SYMBOLIC)