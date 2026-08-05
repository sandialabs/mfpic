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

base_num_elements = 100
number_macroparticles_per_cell = 100
# number_macroparticles_per_cell = 1000000

N2_species = species.Species(mass=4.65e-26, specific_heat_ratio=1.4)
number_density_offset = 1e26
perturbation = 1.0
number_density_height = perturbation * number_density_offset

max_mass_density = (number_density_height + number_density_offset) * N2_species.mass
pressure = 27613
sound_speed = euler.speed_of_sound(N2_species, max_mass_density, pressure)

mach_number = 2.5
velocity = mach_number * sound_speed

final_time = 1 * gaussian_standard_deviation / velocity

def get_input_deck():
    num_elements = base_num_elements

    num_time_steps = 2
    dt = final_time / num_time_steps

    number_macroparticles_per_population = number_macroparticles_per_cell * num_elements

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
      Number of Macroparticles per Species: {number_macroparticles_per_population}
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
  Stride: 1
"""
    return input_deck_contents


def run(mfpic_executable):
    input_deck_contents = get_input_deck()
    yaml = "advecting_gaussian_particles.yaml"
    with open(yaml, "w") as input_deck:
        input_deck.write(input_deck_contents)

    result = subprocess.run([mfpic_executable, "-i", yaml], capture_output=True, text=True)

    log_filename = f"execute.log"
    with open(log_filename, "w") as execute_log:
        execute_log.write(result.stdout)
        execute_log.write(result.stderr)

    result.check_returncode()


def analyze():
    particle_moments = np.genfromtxt('particle_moments_N2.csv', names=True, delimiter=',')
    particle_moments_over_time = dict()
    num_steps = int(np.max(particle_moments['step']))
    for i in range(num_steps + 1):
        mask = particle_moments['step'] == i
        particle_moments_over_time[i] = particle_moments[mask]

    figures_directory = "Figures"
    os.makedirs(figures_directory, exist_ok=True)

    for i in range(num_steps + 1):
        fig, axes = plt.subplots()
        axes.plot(particle_moments_over_time[i]['x'], particle_moments_over_time[i]['number_density'], label="Number Density")
        axes.set_xlabel("x")
        axes.set_ylabel("n")
        fig.savefig(f"{figures_directory}/NumberDensity{i:02}.png")
        plt.close(fig)


if __name__ == "__main__":
    import sys

    if "run" in sys.argv[1:]:
        mfpic_executable = sys.argv[2]
        run(mfpic_executable)
    else:
        analyze()