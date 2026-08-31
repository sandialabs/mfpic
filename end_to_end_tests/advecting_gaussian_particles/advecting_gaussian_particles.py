import sys

sys.path.append("../python")

import euler
import read_mesh_data
import species
import utils

import numpy as np
import matplotlib.pyplot as plt
import subprocess

gaussian_standard_deviation = 0.1
gaussian_center = 4 * gaussian_standard_deviation
domain_length = 9 * gaussian_standard_deviation

base_num_elements = 25
refinement_levels = [4]
# refinement_levels = [16]

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

    mesh_folder_name = format_mesh_folder_name(refinement_level)
    num_macroparticles_per_population = 100 * num_elements

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
  N2:
    Mass: {N2_species.mass}
    Charge: 0
    Specific Heat Ratio: {N2_species.specific_heat_ratio}

Particles:
  Initial Conditions:
    - Species: [N2]
      Number of Macroparticles per Species: {num_macroparticles_per_population}
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
    mesh_folder_name = format_mesh_folder_name(refinement_levels[0])
    timesteps, mesh_data = read_mesh_data.read_mesh_data(mesh_folder_name)
    points = mesh_data[-1]["points"]
    x_points = points[:,0]
    dx = x_points[1] - x_points[0]

    num_cells = base_num_elements * refinement_levels[0]
    x_plot = np.linspace(0, domain_length, 10 * num_cells)
    cell_centers = np.linspace(0.5 * dx, domain_length - 0.5 * dx, num_cells)

    import h5py

    particle_file = h5py.File("particles.h5part")
    particle_data_at_last_timestep = particle_file[f"Step#1"]
    for i, time in enumerate(timesteps):
        particle_data = particle_file[f"Step#{i}"]
        number_of_particles_in_each_cell = np.zeros(num_cells)
        for i_element, weight in zip(particle_data['element'], particle_data['weight']):
            number_of_particles_in_each_cell[i_element] += weight
        mass_density = number_of_particles_in_each_cell * (N2_species.mass / dx)

        fig, axes = plt.subplots()
        axes.plot(cell_centers, mass_density, label="Numerical Solution")
        axes.plot(x_plot, exact_mass_density(x_plot, time), label="Exact Solution")
        axes.legend()
        axes.set_title(f"Mass Density At Time = {time}")
        axes.set_xlabel("x")
        axes.set_ylabel("rho")
        fig.savefig(f"MassDensity{i:02}.png")
        plt.close(fig)


if __name__ == "__main__":
    import sys

    if "run" in sys.argv[1:]:
        run(sys.argv[2])
    elif "plot" in sys.argv[1:]:
        plot()
    else:
        analyze()
