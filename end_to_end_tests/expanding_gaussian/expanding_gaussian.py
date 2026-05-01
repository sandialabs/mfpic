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

final_time = 6.6e-9

max_cfl = 0.3
max_wavespeed = sound_speed
dt, num_time_steps = utils.compute_timestepping_that_satisfies_cfl(max_cfl, dx, max_wavespeed, final_time)

basis_order = 0


def format_mesh_folder_name(time_integrator):
    return f"{time_integrator.replace(" ", "_")}/MeshOutput"


def get_input_deck(time_integrator):
    mesh_folder_name = format_mesh_folder_name(time_integrator)

    input_deck_contents = f"""
Mesh:
  Type: line
  Lengths: [{domain_length}]
  Number of Elements: [{num_elements}]
  Periodic Dimensions: [x]

Time Stepping:
  Number of Time Steps: {num_time_steps}
  Time Step Size: {dt}
  Type: {time_integrator}

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
  Stride: 1
  Mesh Output Folder: {mesh_folder_name}
"""
    return input_deck_contents

def run(mfpic_executable, time_integrator):
    yaml = "expanding_gaussian.yaml"
    input_deck_contents = get_input_deck(time_integrator)
    with open(yaml, "w") as input_deck:
        input_deck.write(input_deck_contents)

    result = subprocess.run([mfpic_executable, "-i", yaml])
    result.check_returncode()

    os.rename("output_lf_0.csv", f"{time_integrator.replace(" ", "_")}/output_lf_0.csv")


def analyze():
    forward_euler_data = np.genfromtxt('Forward_Euler/output_lf_0.csv', names=True)
    verlet_data = np.genfromtxt('Verlet/output_lf_0.csv', names=True)

    time = forward_euler_data['Time']

    field_energy_forward_euler = forward_euler_data['Field_Energy']
    fluid_energy_forward_euler = forward_euler_data['Total_Fluid_Energy']
    total_energy_forward_euler = field_energy_forward_euler + fluid_energy_forward_euler
    fluid_kinetic_energy_forward_euler = forward_euler_data['Total_Fluid_Kinetic_Energy']
    fluid_internal_energy_forward_euler = fluid_energy_forward_euler - fluid_kinetic_energy_forward_euler

    field_energy_verlet = verlet_data['Field_Energy']
    fluid_energy_verlet = verlet_data['Total_Fluid_Energy']
    total_energy_verlet = field_energy_verlet + fluid_energy_verlet
    fluid_kinetic_energy_verlet = verlet_data['Total_Fluid_Kinetic_Energy']
    fluid_internal_energy_verlet = fluid_energy_verlet - fluid_kinetic_energy_verlet

    figures_directory = "Figures"
    os.makedirs(figures_directory, exist_ok=True)

    fig, axes = plt.subplots()
    axes.plot(time, total_energy_forward_euler, label="Forward Euler")
    axes.plot(time, total_energy_verlet, label="Verlet")
    axes.legend()
    axes.set_title(f"Total Energy over Time")
    axes.set_xlabel("time")
    axes.set_ylabel("Energy")
    fig.savefig(f"{figures_directory}/TotalEnergy.png")
    plt.close(fig)

    fig, axes = plt.subplots()
    axes.plot(time, fluid_internal_energy_forward_euler, label="Forward Euler")
    axes.plot(time, fluid_internal_energy_verlet, label="Verlet")
    axes.legend()
    axes.set_title(f"Fluid Internal Energy over Time")
    axes.set_xlabel("time")
    axes.set_ylabel("Internal Energy")
    fig.savefig(f"{figures_directory}/InternalEnergy.png")
    plt.close(fig)

    assert(total_energy_forward_euler[-1] > total_energy_verlet[-1])
    assert(fluid_internal_energy_forward_euler[-1] < fluid_internal_energy_verlet[-1])


def plot(time_integrator):
    mesh_folder_name = format_mesh_folder_name(time_integrator)
    timesteps, mesh_data = read_mesh_data.read_mesh_data(mesh_folder_name)

    points = mesh_data[0]["points"]
    x_points = points[:, 0]

    figures_directory = f"{time_integrator.replace(" ", "_")}/Figures"
    os.makedirs(figures_directory, exist_ok=True)
    for i, time in enumerate(timesteps):
        electron_data = np.transpose(mesh_data[i]["species_0_lf_0"])
        proton_data = np.transpose(mesh_data[i]["species_1_lf_0"])
        electron_mass_density = electron_data[0]
        proton_mass_density = proton_data[0]
        electron_bulk_velocity = euler.get_bulk_velocity_from_conservative_state(electron_data)
        proton_bulk_velocity = euler.get_bulk_velocity_from_conservative_state(proton_data)
        electron_pressure = euler.get_pressure_from_conservative_state(electron_data, electron_species)
        proton_pressure = euler.get_pressure_from_conservative_state(proton_data, proton_species)
        electron_internal_energy_density = euler.get_internal_energy_density_from_conservative_state(electron_data)
        proton_internal_energy_density = euler.get_internal_energy_density_from_conservative_state(proton_data)

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

        fig, axes = plt.subplots()
        axes.plot(x_points, electron_bulk_velocity[0], label="Electrons")
        axes.plot(x_points, proton_bulk_velocity[0], label="Protons")
        axes.legend()
        axes.set_title(f"Bulk Velocity At Time = {time}")
        axes.set_xlabel("x")
        axes.set_ylabel("v")
        fig.savefig(f"{figures_directory}/BulkVelocity{i:02}.png")
        plt.close(fig)

        fig, axes = plt.subplots()
        axes.plot(x_points, electron_pressure, label="Electrons")
        axes.plot(x_points, proton_pressure, label="Protons")
        axes.legend()
        axes.set_title(f"Pressure At Time = {time}")
        axes.set_xlabel("x")
        axes.set_ylabel("P")
        fig.savefig(f"{figures_directory}/Pressure{i:02}.png")
        plt.close(fig)

        fig, axes = plt.subplots()
        axes.plot(x_points, electron_internal_energy_density, label="Electrons")
        axes.plot(x_points, proton_internal_energy_density, label="Protons")
        axes.legend()
        axes.set_title(f"Internal Energy Density At Time = {time}")
        axes.set_xlabel("x")
        axes.set_ylabel("P")
        fig.savefig(f"{figures_directory}/InternalEnergyDensity{i:02}.png")
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

    txt_data = np.genfromtxt(f"{time_integrator.replace(" ", "_")}/output_lf_0.csv", names=True)
    time = txt_data['Time']
    field_energy = txt_data['Field_Energy']
    fluid_energy = txt_data['Total_Fluid_Energy']
    total_energy = field_energy + fluid_energy
    fluid_kinetic_energy = txt_data['Total_Fluid_Kinetic_Energy']
    fluid_internal_energy = fluid_energy - fluid_kinetic_energy

    fig, axes = plt.subplots()
    axes.plot(time, total_energy, label="Total Energy")
    axes.plot(time, field_energy, label="Field Energy")
    axes.plot(time, fluid_energy, label="Fluid Energy")
    axes.legend()
    axes.set_title(f"Energies Over Time")
    axes.set_xlabel("time")
    axes.set_ylabel("Energy")
    fig.savefig(f"{figures_directory}/Energies.png")
    plt.close(fig)

    fig, axes = plt.subplots()
    axes.plot(time, fluid_energy, label="Total Fluid Energy")
    axes.plot(time, fluid_kinetic_energy, label="Fluid Kinetic Energy")
    axes.plot(time, fluid_internal_energy, label="Fluid Internal Energy")
    axes.legend()
    axes.set_title(f"Fluid Energies Over Time")
    axes.set_xlabel("time")
    axes.set_ylabel("Energy")
    fig.savefig(f"{figures_directory}/FluidEnergies.png")
    plt.close(fig)

if __name__ == "__main__":
    import sys

    time_integrators = ["Forward Euler", "Verlet"]
    if "run" in sys.argv[1:]:
        for time_integrator in time_integrators:
            run(sys.argv[2], time_integrator)
    if "plot" in sys.argv[1:]:
        for time_integrator in time_integrators:
            plot(time_integrator)
    else:
        analyze()
