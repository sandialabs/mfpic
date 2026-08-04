import sys
sys.path.append("../python")

import plasma_parameters
import species

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

final_time = plasma_period
num_time_steps = 10
dt = final_time / num_time_steps

domain_length = debye_length
num_elements = 10
dx = domain_length / num_elements

number_macroparticles_per_cell = 1000
number_macroparticles_per_population = number_macroparticles_per_cell * num_elements

def get_input_deck():
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
"""
    return input_deck_contents


def run(mfpic_executable):
    input_deck_contents = get_input_deck()
    yaml = "spatially_constant_maxwellian.yaml"
    with open(yaml, "w") as input_deck:
        input_deck.write(input_deck_contents)

    result = subprocess.run([mfpic_executable, "-i", yaml], capture_output=True, text=True)

    log_filename = f"execute.log"
    with open(log_filename, "w") as execute_log:
        execute_log.write(result.stdout)
        execute_log.write(result.stderr)

    result.check_returncode()
    pass

def analyze():
    particle_moments = np.genfromtxt('particle_moments_electron.csv', names=True, delimiter=',')
    particle_moments_over_time = dict()
    num_steps = int(np.max(particle_moments['step']))
    for i in range(num_steps + 1):
        mask = particle_moments['step'] == i
        particle_moments_over_time[i] = particle_moments[mask]

    figures_directory = "Figures"
    os.makedirs(figures_directory, exist_ok=True)

    for i in range(num_steps + 1):
        fig, axes = plt.subplots()
        x_points = particle_moments_over_time[i]['x']
        axes.plot(x_points, particle_moments_over_time[i]['number_density'], label="MFPIC")
        axes.plot(x_points, np.ones(x_points.shape) * number_density, label="EMPIRE")
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