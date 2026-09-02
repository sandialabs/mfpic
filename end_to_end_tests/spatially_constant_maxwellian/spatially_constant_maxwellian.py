import sys
sys.path.append("../python")

from compute_particle_distribution_function import *
import maxwellian
import plasma_parameters
from plot_particle_distribution_function import *
import species

import h5py
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.constants import electron_mass, proton_mass, elementary_charge
import subprocess

number_density = 1e16
temperature = 300

electron_species = species.Species(mass=electron_mass, charge=-elementary_charge)
proton_species = species.Species(mass=proton_mass, charge=elementary_charge)

debye_length = plasma_parameters.debye_length(electron_species, number_density, temperature)
plasma_frequency = plasma_parameters.plasma_frequency(electron_species, number_density)
plasma_period = plasma_parameters.plasma_period(electron_species, number_density)

num_plasma_periods = 10
final_time = num_plasma_periods * plasma_period
num_time_steps = 20 * num_plasma_periods
dt = final_time / num_time_steps

domain_length = 100 * debye_length
debye_lengths_per_cell = 10
dx = debye_lengths_per_cell * debye_length
num_cells = int(domain_length / dx)

R = 0.5 * dx
R_over_debye_length = R / debye_length


def exact_particle_distribution_function(x, v, t):
    bulk_velocity = 0
    return number_density * maxwellian.probability_density_function_1d(v, electron_species, bulk_velocity, temperature)


def exact_number_density(x, t):
    return number_density


def get_output_directory(number_macroparticles_per_cell, debye_lengths_per_cell):
    return f"{number_macroparticles_per_cell}ParticlesPerCell{debye_lengths_per_cell}DebyeLengthsPerCell"


def get_input_deck(number_macroparticles_per_cell):
    output_directory = get_output_directory(number_macroparticles_per_cell, debye_lengths_per_cell)
    number_macroparticles_per_population = number_macroparticles_per_cell * num_cells
    input_deck_contents = f"""
Mesh:
  Type: line
  Lengths: [{domain_length}]
  Number of Elements: [{num_cells}]
  Periodic Dimensions: [x]

Time Stepping:
  Number of Time Steps: {num_time_steps}
  Time Step Size: {dt}

Species:
  electron:
    Mass: {electron_species.mass}
    Charge: {electron_species.charge}
    Specific Heat Ratio: {electron_species.specific_heat_ratio}
  proton:
    Mass: {proton_species.mass}
    Charge: {proton_species.charge}
    Specific Heat Ratio: {proton_species.specific_heat_ratio}

Particles:
  Initial Conditions:
    - Species: [electron, proton]
      Number of Macroparticles per Species: {number_macroparticles_per_population}
      Constant:
        Number Density: {number_density}
        Temperature: {temperature}

Output:
  Stride: 1
  Mesh Output Folder: {output_directory}/MeshOutput
  Particle Dump Filename: particles.h5part
"""
    return input_deck_contents


def compute_particle_weight(number_macroparticles_per_cell):
    particle_weight = (number_density * np.power(debye_length, 3) / number_macroparticles_per_cell) * (dx / debye_length)
    return particle_weight


def run(mfpic_executable, number_macroparticles_per_cell):
    output_directory = get_output_directory(number_macroparticles_per_cell, debye_lengths_per_cell)
    input_deck_contents = get_input_deck(number_macroparticles_per_cell)
    yaml = "spatially_constant_maxwellian.yaml"
    with open(yaml, "w") as input_deck:
        input_deck.write(input_deck_contents)

    cell_edges = np.linspace(0, domain_length, num_cells + 1)
    thermal_speed = maxwellian.thermal_speed(electron_species, temperature)
    min_velocity = -30 * thermal_speed
    max_velocity = 30 * thermal_speed
    velocity_bins = np.linspace(min_velocity, max_velocity, 200)
    times = np.linspace(0, final_time, num_time_steps + 1)
    empty_pdf_data = generate_empty_particle_distribution_function_data(cell_edges, velocity_bins, times)
    with h5py.File("electron_particle_distribution_function_data.h5", 'w') as hf:
        write_dict_to_hdf5_file(hf, empty_pdf_data)

    num_runs = 500
    for i_run in range(num_runs):
        result = subprocess.run([mfpic_executable, "-i", yaml], capture_output=True, text=True)

        log_filename = f"{output_directory}/execute.log"
        with open(log_filename, "w") as execute_log:
            execute_log.write(result.stdout)
            execute_log.write(result.stderr)

        result.check_returncode()

        with h5py.File("particles.h5part") as particle_data:
            particle_moment_data = np.genfromtxt('particle_moments_electron.csv', names=True, delimiter=',')
            particle_moment_data = arrange_particle_moment_data_by_timestep(particle_moment_data)

            with h5py.File("electron_particle_distribution_function_data.h5", 'r+') as particle_distribution_function_data:
                update_particle_distribution_function_data(
                    particle_distribution_function_data,
                    particle_data,
                    particle_moment_data,
                    b'electron')

        # files_to_copy = [yaml, "particle_moments_electron.csv", "particle_moments_proton.csv", "output.csv", "particles.h5part"]
        # for file in files_to_copy:
        #     os.rename(file, f"{output_directory}/{file}")


def analyze_single_run(number_macroparticles_per_cell):
    output_directory = get_output_directory(number_macroparticles_per_cell, debye_lengths_per_cell)
    error_over_time = compute_error_over_time(number_macroparticles_per_cell)
    particle_weight = compute_particle_weight(number_macroparticles_per_cell)

    figures_directory = f"{output_directory}/Figures"
    os.makedirs(figures_directory, exist_ok=True)

    times = np.linspace(0, final_time / plasma_period, num_time_steps+1)
    fig, axes = plt.subplots()
    axes.plot(times, error_over_time)
    axes.set_xlabel("Plasma Periods")
    axes.set_ylabel("Error in f")
    axes.set_title(f"R/debye_length = {R_over_debye_length}, W = {particle_weight:.1e}")
    fig.savefig(f"{figures_directory}/fError.png")
    plt.close(fig)


def analyze(particles_per_cell_array):
    with h5py.File("electron_particle_distribution_function_data.h5", 'r') as particle_distribution_function:
        plot_particle_distribution_function(particle_distribution_function, exact_particle_distribution_function)

        thermal_speed = maxwellian.thermal_speed(electron_species, temperature)
        v_max = 10 * thermal_speed
        plot_particle_distribution_function(particle_distribution_function, exact_particle_distribution_function, v_max=v_max)

        plot_number_density(particle_distribution_function, exact_number_density)

        with h5py.File("particles.h5part", 'r') as particle_data:
            plot_particle_distribution_function_vs_stored_f(particle_data, particle_distribution_function, b'electron', exact_particle_distribution_function)

            times = np.linspace(0, final_time, num_time_steps + 1)
            error_in_stored_f_vs_exact_f_over_time = compute_error_in_stored_f_vs_exact_f(particle_data, b'electron', exact_particle_distribution_function, times)
            error_in_stored_f_vs_histogram_f_over_time = compute_error_in_stored_f_vs_histogram_f(particle_data, b'electron', particle_distribution_function)

            fig, axes = plt.subplots()
            axes.plot(times, error_in_stored_f_vs_exact_f_over_time, label="particle_f - exact_f")
            axes.plot(times, error_in_stored_f_vs_histogram_f_over_time, label="particle_f - histogram_f")
            axes.legend()
            axes.set_xlabel("time")
            axes.set_ylabel("np.sum(diff) / num_macroparticles")
            fig.savefig(f"Figures/ErrorOverTime.png")
            plt.close(fig)


        # loop over particles per cell
        # loop over dx size

    # error over time between exact particle distribution function and ensemble_averaged
    # difference between f stored on particles and ensemble averaged

    # figures_directory = f"Figures"
    # os.makedirs(figures_directory, exist_ok=True)

    # errors_over_time_dict = dict()
    # for number_macroparticles_per_cell in particles_per_cell_array:
    #     errors_over_time_dict[number_macroparticles_per_cell] = compute_error_over_time(number_macroparticles_per_cell)

    # fig, axes = plt.subplots()
    # times = np.linspace(0, final_time / plasma_period, num_time_steps+1)
    # for number_macroparticles_per_cell in particles_per_cell_array:
    #     particle_weight = compute_particle_weight(number_macroparticles_per_cell)
    #     axes.plot(times, errors_over_time_dict[number_macroparticles_per_cell], label=f"W = {particle_weight:.1e}")
    # axes.legend()
    # axes.set_xlabel("Plasma Periods")
    # axes.set_ylabel("Error in f")
    # axes.set_title(f"R/debye_length = {R_over_debye_length}")
    # fig.savefig(f"{figures_directory}/fError{debye_lengths_per_cell}DebyeLengthsPerCell.png")
    # plt.close(fig)


if __name__ == "__main__":
    import sys

    # particles_per_cell_array = [10, 100, 1000]
    particles_per_cell_array = [1000]
    if "run" in sys.argv[1:]:
        mfpic_executable = sys.argv[2]
        for number_macroparticles_per_cell in particles_per_cell_array:
            run(mfpic_executable, number_macroparticles_per_cell)
    else:
        # for number_macroparticles_per_cell in particles_per_cell_array:
        #     analyze_single_run(number_macroparticles_per_cell)

        analyze(particles_per_cell_array)
