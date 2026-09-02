import sys

sys.path.append("../python")

from compute_particle_distribution_function import *
from plot_particle_distribution_function import *
import euler
import maxwellian
import read_mesh_data
import species
import utils
import verification

import h5py
import matplotlib.pyplot as plt
import numpy as np
import os
import re
from scipy.constants import Boltzmann as boltzmann_constant
import scipy.integrate
import subprocess

gaussian_standard_deviation = 0.1
gaussian_center = 4 * gaussian_standard_deviation
domain_length = 9 * gaussian_standard_deviation

base_num_cells = 100
num_cells = base_num_cells
number_macroparticles_per_cell = 10000
# number_macroparticles_per_cell = 1000000

N2_species = species.Species(mass=4.65e-26, specific_heat_ratio=1.4)
number_density_offset = 1e26
perturbation = 1.0
number_density_height = perturbation * number_density_offset

max_mass_density = (number_density_height + number_density_offset) * N2_species.mass
pressure = 27613
sound_speed = euler.speed_of_sound(N2_species, max_mass_density, pressure)

mach_number = 2.5
bulk_velocity = mach_number * sound_speed

final_time = 9 * gaussian_standard_deviation / bulk_velocity
num_time_steps = 2
dt = final_time / num_time_steps

number_macroparticles_per_population = number_macroparticles_per_cell * num_cells

def initial_number_density(x):
    exponent = -0.5 * np.power((x - gaussian_center) / gaussian_standard_deviation, 2)
    return number_density_height * np.exp(exponent) + number_density_offset


def initial_f_v(x, v):
    temperature = pressure / (initial_number_density(x) * boltzmann_constant)
    return maxwellian.probability_density_function_1d(v, N2_species, bulk_velocity, temperature)


def initial_particle_distribution_function(x, v):
    return initial_number_density(x) * initial_f_v(x, v)


def exact_particle_distribution_function(x, v, t):
    initial_x = np.mod(x - v * t, domain_length)
    initial_v = v
    return initial_particle_distribution_function(initial_x, initial_v)


def exact_number_density(x, t):
    f = lambda v : exact_particle_distribution_function(x, v, t)
    integral, _ = scipy.integrate.quad(f, -np.inf, np.inf, epsrel=1e-8, limit=100)
    return integral


def get_input_deck():

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
  N2:
    Mass: {N2_species.mass}
    Charge: 0
    Specific Heat Ratio: {N2_species.specific_heat_ratio}

Particles:
  Initial Conditions:
    - Species: [N2]
      Number of Macroparticles per Species: {number_macroparticles_per_population}
      Gaussian:
        Center: [{gaussian_center}]
        Standard Deviation: {gaussian_standard_deviation}
        Offsets:
          Number Density: {number_density_offset}
          Pressure: {pressure}
          Bulk Velocity: [{bulk_velocity}, 0., 0.]
        Heights:
          Number Density: {number_density_height}
          Pressure: 0.

Output:
  Stride: 1
"""
    return input_deck_contents


def run(mfpic_executable):
    input_deck_contents = get_input_deck()
    yaml = "advecting_gaussian_particles.yaml"
    with open(yaml, "w") as input_deck:
        input_deck.write(input_deck_contents)

    cell_edges = np.linspace(0, domain_length, num_cells + 1)
    temperature = pressure / (number_density_offset * boltzmann_constant)
    thermal_speed = maxwellian.thermal_speed(N2_species, temperature)
    min_velocity = bulk_velocity - 5.6 * thermal_speed
    max_velocity = bulk_velocity + 5.6 * thermal_speed
    velocity_bins = np.linspace(min_velocity, max_velocity, 100)
    times = np.linspace(0, final_time, num_time_steps + 1)
    empty_pdf_data = generate_empty_particle_distribution_function_data(cell_edges, velocity_bins, times)
    with h5py.File("particle_distribution_function_data.h5", 'w') as hf:
        write_dict_to_hdf5_file(hf, empty_pdf_data)

    num_runs = 2
    for i_run in range(num_runs):
        print(f"i_run = {i_run}")
        result = subprocess.run([mfpic_executable, "-i", yaml], capture_output=True, text=True)

        log_filename = f"execute.log"
        with open(log_filename, "w") as execute_log:
            execute_log.write(result.stdout)
            execute_log.write(result.stderr)

        result.check_returncode()

        with h5py.File("particles.h5part") as particle_data:
            particle_moment_data = np.genfromtxt('particle_moments_N2.csv', names=True, delimiter=',')
            particle_moment_data = arrange_particle_moment_data_by_timestep(particle_moment_data)
            with h5py.File("particle_distribution_function_data.h5", 'r+') as particle_distribution_function_data:
                update_particle_distribution_function_data(
                    particle_distribution_function_data,
                    particle_data,
                    particle_moment_data,
                    b'N2')


def plot_particle_stored_f_vs_exact_f(particle_data):
    min_v = np.min([np.min(particle_data_i['vx']) for particle_data_i in particle_data.values()])
    max_v = np.max([np.max(particle_data_i['vx']) for particle_data_i in particle_data.values()])
    v_plot = np.linspace(min_v, max_v, 200)

    for i, (timestep_key, particle_data_i) in enumerate(particle_data.items()):
        figures_directory = f"Figures/{timestep_key}/ParticleStoredfVsExactf"
        os.makedirs(figures_directory, exist_ok=True)
        time = i * dt

        exact_f = exact_particle_distribution_function(particle_data_i['x'], particle_data_i['vx'], time)
        errors = np.abs(particle_data_i['particle_distribution_function_value'] - exact_f) / exact_f
        print(np.linalg.norm(errors))
        fig, axes = plt.subplots()
        axes.plot(np.array(particle_data_i['vx']), errors, 'o')
        axes.set_xlabel("v")
        axes.set_ylabel("f Error")
        fig.savefig(f"{figures_directory}/Errors.png")
        plt.close(fig)

        for i_cell in range(num_cells):
            mask = np.array(particle_data_i['element']) == i_cell
            velocities = particle_data_i['vx'][mask]
            positions = particle_data_i['x'][mask]
            f = particle_data_i['particle_distribution_function_value'][mask]

            fig, axes = plt.subplots()
            axes.plot(velocities, f, 'o', label="Particles")
            for x in positions[:4]:
                axes.plot(v_plot, exact_particle_distribution_function(x, v_plot, time), label=f"x = {x:.3e}")
            axes.legend()
            axes.set_xlabel("v")
            axes.set_ylabel("f")
            fig.savefig(f"{figures_directory}/Cell{i_cell:03d}.png")
            plt.close(fig)


def analyze():
    with h5py.File("particles.h5part") as particle_data:
        with h5py.File("particle_distribution_function_data.h5", 'r') as particle_distribution_function_data:
            plot_number_density(particle_data, particle_distribution_function_data)
            plot_particle_distribution_function(particle_distribution_function_data, exact_particle_distribution_function)
            plot_particle_stored_f_vs_exact_f(particle_data)


if __name__ == "__main__":
    import sys

    if "run" in sys.argv[1:]:
        mfpic_executable = sys.argv[2]
        run(mfpic_executable)
    else:
        analyze()