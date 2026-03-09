import sys

sys.path.append("../python")

import euler
import read_mesh_data
import species
import utils
import verification

import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess

gaussian_standard_deviation = 0.1
gaussian_center = 4 * gaussian_standard_deviation
domain_length = 9 * gaussian_standard_deviation

base_num_elements = 25
refinement_levels = [16, 32, 64]

N2_species = species.Species(mass=4.65e-26, specific_heat_ratio=1.4)
number_density_offset = 1e26
perturbation = 0.001
number_density_height = perturbation * number_density_offset

max_mass_density = (number_density_height + number_density_offset) * N2_species.mass
pressure = 27613
sound_speed = euler.speed_of_sound(N2_species, max_mass_density, pressure)

mach_number = 2.5
velocity = mach_number * sound_speed

max_cfl = 0.95
max_wavespeed = velocity + sound_speed

final_time = 10 * gaussian_standard_deviation / velocity

basis_order = 0

def exact_mass_density(x, t):
    sigma = gaussian_standard_deviation
    c = gaussian_center
    mask = x < np.fmod(velocity * t, domain_length)
    number_of_periods = int(velocity * t / domain_length)
    shift = (mask + number_of_periods) * domain_length
    exponential = np.exp(-0.5 * np.power(shift + x - velocity * t - c, 2) / np.power(sigma, 2))
    number_density = number_density_height * exponential + number_density_offset
    return N2_species.mass * number_density


def format_mesh_folder_name(refinement_level):
    return f"MeshOutput{refinement_level:02}"


def get_input_deck(refinement_level):
    num_elements = base_num_elements * refinement_level
    dx = domain_length / num_elements
    dt, num_time_steps = utils.compute_timestepping_that_satisfies_cfl(max_cfl, dx, max_wavespeed, final_time)

    # num_time_steps = 1

    mesh_folder_name = format_mesh_folder_name(refinement_level)

    input_deck_contents = f"""
Mesh:
  Type: line
  Lengths: [{domain_length}]
  Number of Elements: [{num_elements}]
  Periodic Dimensions: [x]

Time Stepping:
  Number of Time Steps: {num_time_steps}
  Time Step Size: {dt}

Species:
  neutral_electron:
    Mass: {N2_species.mass}
    Charge: 0
    Specific Heat Ratio: {N2_species.specific_heat_ratio}

Euler Fluids:
  Basis Order: {basis_order}
  Initial Conditions:
    - Species: [neutral_electron]
      Gaussian:
        Center: [{gaussian_center}]
        Standard Deviation: {gaussian_standard_deviation}
        Offsets:
          Number Density: {number_density_offset}
          Pressure: {pressure}
          Bulk Velocity: [{velocity}, 0., 0.]
        Heights:
          Number Density: {number_density_height}
          Pressure: 0.

Output:
  Stride: {num_time_steps}
  Mesh Output Folder: {mesh_folder_name}
"""
    return input_deck_contents


def run(mfpic_executable):
    for refinement_level in refinement_levels:
        input_deck_contents = get_input_deck(refinement_level)
        yaml = "advecting_gaussian.yaml"
        with open(yaml, "w") as input_deck:
            input_deck.write(input_deck_contents)

        print(mfpic_executable, "-i", yaml)
        result = subprocess.run([mfpic_executable, "-i", yaml])
        result.check_returncode()


def compute_error(data, points, exact_solution):
    numerical_solution = verification.create_1D_interpolater(data, points)
    exact_solution_at_final_time = lambda x: exact_solution(x, final_time)
    error = verification.compute_L1_relative_error_1D(numerical_solution, exact_solution_at_final_time, points)
    return error


def analyze():
    errors = []
    h_list = []
    for refinement_level in refinement_levels:
        print(f"refinement_level = {refinement_level}")
        mesh_folder_name = format_mesh_folder_name(refinement_level)
        print("read mesh data")
        _, mesh_data = read_mesh_data.read_mesh_data(mesh_folder_name)
        points = mesh_data[-1]["points"]
        x_points = points[:,0]
        dx = x_points[1] - x_points[0]
        print(f"dx = {dx}")
        h_list.append(dx)

        fluid_data = np.transpose(mesh_data[-1]["species_0"])
        print("compute error")
        error = compute_error(fluid_data[0], x_points, exact_mass_density)
        print(f"error = {error}")
        errors.append(error)

    rates = verification.compute_convergence_rates(errors, h_list)
    print(f"rates = {rates}")

    verification.plot_errors_and_expected_convergence_rate(h_list, errors, 1.0)

    expected_convergence_rate = basis_order + 1
    tolerance = 0.1
    assert(np.all(rates > expected_convergence_rate - tolerance))


def plot():
    for refinement_level in refinement_levels:
        mesh_folder_name = format_mesh_folder_name(refinement_level)
        timesteps, mesh_data = read_mesh_data.read_mesh_data(mesh_folder_name)

        points = mesh_data[0]["points"]
        x_points = points[:, 0]
        num_cells = int(0.5 * x_points.shape[0])
        x_plot = np.linspace(0, domain_length, 10 * num_cells)

        figures_directory = f"Figures{refinement_level:02}"
        os.makedirs(figures_directory, exist_ok=True)
        for i, time in enumerate(timesteps):
            fluid_data = np.transpose(mesh_data[i]["species_0"])
            mass_density_data = fluid_data[0]

            fig, axes = plt.subplots()
            axes.plot(x_points, mass_density_data, label="Numerical Solution")
            axes.plot(x_plot, exact_mass_density(x_plot, time), label="Exact Solution")
            axes.legend()
            axes.set_title(f"Mass Density At Time = {time}")
            axes.set_xlabel("x")
            axes.set_ylabel("rho")
            fig.savefig(f"{figures_directory}/MassDensity{i:02}.png")
            plt.close(fig)


if __name__ == "__main__":
    import sys

    if "run" in sys.argv[1:]:
        run(sys.argv[2])
    elif "plot" in sys.argv[1:]:
        plot()
    else:
        analyze()
