import sys

sys.path.append("../python")

import euler
import read_mesh_data
import species
import utils

import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.constants import Boltzmann, electron_mass, proton_mass, elementary_charge
import subprocess

domain_length = 0.01
gaussian_standard_deviation = 0.0005
gaussian_center = domain_length / 2

num_elements = 300
dx = domain_length / num_elements

proton_species = species.Species(charge=elementary_charge, mass=proton_mass)
electron_species = species.Species(charge=-elementary_charge, mass=electron_mass)

number_density_peak = 1e16
number_density_offset = 1. / 8. * number_density_peak
number_density_height = number_density_peak - number_density_offset

temperature_peak = 1000
temperature_offset = 0.1 * temperature_peak

pressure_peak = number_density_peak * Boltzmann * temperature_peak
pressure_offset = number_density_offset * Boltzmann * temperature_offset
pressure_height = pressure_peak - pressure_offset

electron_mass_density_peak = number_density_peak * electron_species.mass
sound_speed = euler.speed_of_sound(electron_species, electron_mass_density_peak, pressure_peak)

final_time = 1e-7

max_cfl = 0.5
max_wavespeed = sound_speed
dt, num_time_steps = utils.compute_timestepping_that_satisfies_cfl(max_cfl, dx, max_wavespeed, final_time)
# num_time_steps = 5

basis_order = 0

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
  proton:
    Mass: {proton_species.mass}
    Charge: {proton_species.charge}

Euler Fluids:
  Basis Order: {basis_order}
  Initial Conditions:
    - Species: [electron, proton]
      Gaussian:
        Center: [{gaussian_center}]
        Standard Deviation: {gaussian_standard_deviation}
        Offsets:
          Number Density: {number_density_offset}
          Pressure: {pressure_offset}
        Heights:
          Number Density: {number_density_height}
          Pressure: {pressure_height}

Output:
  Stride: 50
"""


def run(mfpic_executable):
    yaml = "expanding_gaussian.yaml"
    with open(yaml, "w") as input_deck:
        input_deck.write(input_deck_contents)

    result = subprocess.run([mfpic_executable, "-i", yaml])
    result.check_returncode()


def analyze():
    pass


def plot():
    mesh_folder_name = "MeshOutput"
    timesteps, mesh_data = read_mesh_data.read_mesh_data(mesh_folder_name)

    points = mesh_data[0]["points"]
    x_points = points[:, 0]

    figures_directory = f"Figures"
    os.makedirs(figures_directory, exist_ok=True)
    for i, time in enumerate(timesteps):
        electron_data = np.transpose(mesh_data[i]["species_0_lf_0"])
        proton_data = np.transpose(mesh_data[i]["species_1_lf_0"])
        electron_mass_density = electron_data[0]
        proton_mass_density = proton_data[0]

        electron_number_density = electron_mass_density / electron_species.mass
        proton_number_density = proton_mass_density / proton_species.mass

        fig, axes = plt.subplots()
        axes.plot(x_points, electron_number_density, label="Electrons")
        axes.plot(x_points, proton_number_density, label="Protons")
        axes.legend()
        axes.set_title(f"Number Density At Time = {time}")
        axes.set_xlabel("x")
        axes.set_ylabel("n")
        fig.savefig(f"{figures_directory}/NumberDensity{i:02}.png")
        plt.close(fig)

        charge_density = (proton_number_density - electron_number_density) * elementary_charge
        fig, axes = plt.subplots()
        axes.plot(x_points, charge_density, label="Charge Density")
        axes.legend()
        axes.set_title(f"Charge Density At Time = {time}")
        axes.set_xlabel("x")
        axes.set_ylabel("rho")
        fig.savefig(f"{figures_directory}/ChargeDensity{i:02}.png")
        plt.close(fig)

        potential = mesh_data[i]["electrostatic_potential_lf_0"]
        fig, axes = plt.subplots()
        axes.plot(x_points, potential, label="Potential")
        axes.legend()
        axes.set_title(f"Potential At Time = {time}")
        axes.set_xlabel("x")
        axes.set_ylabel("Phi")
        fig.savefig(f"{figures_directory}/Potential{i:02}.png")
        plt.close(fig)

        e_field = mesh_data[i]["E_0_lf_0"]
        fig, axes = plt.subplots()
        axes.plot(x_points, e_field, label="E")
        axes.legend()
        axes.set_title(f"E At Time = {time}")
        axes.set_xlabel("x")
        axes.set_ylabel("E")
        fig.savefig(f"{figures_directory}/E{i:02}.png")
        plt.close(fig)



if __name__ == "__main__":
    import sys

    if "run" in sys.argv[1:]:
        run(sys.argv[2])
    if "plot" in sys.argv[1:]:
        plot()
    else:
        analyze()
