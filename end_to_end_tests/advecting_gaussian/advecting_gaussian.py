import sys

sys.path.append("../python")

import euler
import read_mesh_data
import species
import utils

import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess

domain_length = 2.0
base_num_elements = 500
refinement_levels = [2]

gaussian_center = 0.25 * domain_length
gaussian_standard_deviation = 0.1

number_density_offset = 1e26
perturbation = 0.01
number_density_height = perturbation * number_density_offset

pressure = 1000
velocity = 3.0
max_cfl = 0.8

final_time = 0.5 * domain_length / velocity

neutral_species = species.Species(mass=4.65e-26, specific_heat_ratio=1.4)
max_mass_density = (number_density_height + number_density_offset) * neutral_species.mass
sound_speed = euler.speed_of_sound(neutral_species, max_mass_density, pressure)
max_wavespeed = velocity + sound_speed


def exact_mass_density(x, t):
    sigma = gaussian_standard_deviation
    c = gaussian_center
    exponential = np.exp(-0.5 * np.power(x - c - velocity * t, 2) / np.power(sigma, 2))
    number_density = number_density_height * exponential + number_density_offset
    return neutral_species.mass * number_density


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
    Mass: {neutral_species.mass}
    Charge: 0
    Specific Heat Ratio: {neutral_species.specific_heat_ratio}

Euler Fluids:
  Basis Order: 0
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


def analyze():
    pass


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
