import sys

sys.path.append("../python")

import euler
import read_mesh_data
import species
import utils

import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.constants import Boltzmann, electron_mass, elementary_charge, epsilon_0
import subprocess

domain_length = 0.02
gaussian_standard_deviation = 0.0005
gaussian_center = 0.015

num_elements = 600
dx = domain_length / num_elements

proton_mass = 1836 * electron_mass
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

final_time = 2e-6

max_cfl = 0.3
max_wavespeed = sound_speed
dt, num_time_steps = utils.compute_timestepping_that_satisfies_cfl(max_cfl, dx, max_wavespeed, final_time)

# dt = 6.291840740927167e-11
# num_time_steps = 63409

basis_order = 0

plasma_frequency = np.sqrt(number_density_peak * np.power(electron_species.charge, 2) / (electron_mass * epsilon_0))
plasma_period = 2.0 * np.pi / plasma_frequency

print(f"plasma_period = {plasma_period}")
print(f"dt = {dt}")


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
  Stride: 200
  Mesh Output Folder: {mesh_folder_name}
"""
    return input_deck_contents

def run(mfpic_executable, time_integrator):
    assert(10 * dt < plasma_period)

    yaml = "expanding_gaussian.yaml"
    input_deck_contents = get_input_deck(time_integrator)
    with open(yaml, "w") as input_deck:
        input_deck.write(input_deck_contents)

    result = subprocess.run([mfpic_executable, "-i", yaml], capture_output=True, text=True)

    log_filename = "execute.log"
    with open(log_filename, "w") as execute_log:
        execute_log.write(result.stdout)
        execute_log.write(result.stderr)

    os.rename(log_filename, f"{time_integrator.replace(" ", "_")}/{log_filename}")
    os.rename(yaml, f"{time_integrator.replace(" ", "_")}/{yaml}")
    os.rename("output_lf_0.csv", f"{time_integrator.replace(" ", "_")}/output_lf_0.csv")
    os.rename("solver.log", f"{time_integrator.replace(" ", "_")}/solver.log")

    result.check_returncode()


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

    fig, ax = plt.subplots()
    ax.plot(time, total_energy_forward_euler, label="Forward Euler")
    ax.plot(time, total_energy_verlet, label="Verlet")
    ax.legend()
    ax.set_title(f"Total Energy over Time")
    ax.set_xlabel("time")
    ax.set_ylabel("Energy")
    fig.savefig(f"{figures_directory}/TotalEnergy.png")
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.plot(time, fluid_internal_energy_forward_euler, label="Forward Euler")
    ax.plot(time, fluid_internal_energy_verlet, label="Verlet")
    ax.legend()
    ax.set_title(f"Fluid Internal Energy over Time")
    ax.set_xlabel("time")
    ax.set_ylabel("Internal Energy")
    fig.savefig(f"{figures_directory}/InternalEnergy.png")
    plt.close(fig)

    # assert(total_energy_forward_euler[-1] > total_energy_verlet[-1])
    # assert(fluid_internal_energy_forward_euler[-1] < fluid_internal_energy_verlet[-1])


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
        charge_density = (proton_number_density - electron_number_density) * elementary_charge

        potential = mesh_data[i]["electrostatic_potential_lf_0"]
        e_field = mesh_data[i]["E_0_lf_0"]

        fig, ax = plt.subplots()
        ax.plot(x_points, electron_number_density, label="Electrons")
        ax.plot(x_points, proton_number_density, label="Protons")
        ax.legend(loc='upper right')
        ax.set_title(f"Number Density At Time = {time}")
        ax.set_xlabel("x")
        ax.set_ylabel("n")
        ax.set_ylim(0, 1.01 * number_density_peak)
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/NumberDensity{i:03}.png")
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.plot(x_points, electron_bulk_velocity[0], label="Electrons")
        ax.plot(x_points, proton_bulk_velocity[0], label="Protons")
        ax.legend(loc='upper right')
        ax.set_title(f"Bulk Velocity At Time = {time}")
        ax.set_xlabel("x")
        ax.set_ylabel("v")
        ax.set_ylim(-10000, 10000)
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/BulkVelocity{i:03}.png")
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.plot(x_points, electron_pressure, label="Electrons")
        ax.plot(x_points, proton_pressure, label="Protons")
        ax.legend(loc='upper right')
        ax.set_title(f"Pressure At Time = {time}")
        ax.set_xlabel("x")
        ax.set_ylabel("P")
        ax.set_ylim(0, 2.5e-4)
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/Pressure{i:03}.png")
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.plot(x_points, electron_internal_energy_density, label="Electrons")
        ax.plot(x_points, proton_internal_energy_density, label="Protons")
        ax.legend(loc='upper right')
        ax.set_title(f"Internal Energy Density At Time = {time}")
        ax.set_xlabel("x")
        ax.set_ylabel("e")
        ax.set_ylim(0, 3.5e-4)
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/InternalEnergyDensity{i:03}.png")
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.plot(x_points, charge_density, label="Charge Density")
        # ax.legend()
        ax.set_title(f"Charge Density At Time = {time}")
        ax.set_xlabel("x")
        ax.set_ylabel("rho")
        ax.set_ylim([-1.5e-4, 2.5e-4])
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/ChargeDensity{i:03}.png")
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.plot(x_points, potential, label="Potential")
        # ax.legend()
        ax.set_title(f"Potential At Time = {time}")
        ax.set_xlabel("x")
        ax.set_ylabel("Phi")
        ax.set_ylim(-0.35, 0.25)
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/Potential{i:03}.png")
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.plot(x_points, e_field, label="E")
        # ax.legend()
        ax.set_title(f"Electric Field At Time = {time}")
        ax.set_xlabel("x")
        ax.set_ylabel("E")
        ax.set_ylim(-1300, 1300)
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/ElectricField{i:03}.png")
        plt.close(fig)

        fig, axs = plt.subplots(2, 2, figsize=(16, 9))
        axs[0, 0].plot(x_points, electron_number_density, label="Electrons")
        axs[0, 0].plot(x_points, proton_number_density, label="Protons")
        axs[0, 0].legend(loc='upper right')
        axs[0, 0].set_title(f"Number Density")
        axs[0, 0].set_xticklabels([])
        axs[0, 0].set_ylabel("n")
        axs[0, 0].set_ylim(0, 1.01 * number_density_peak)
        axs[0, 1].plot(x_points, charge_density, label="Charge Density")
        axs[0, 1].set_title(f"Charge Density")
        axs[0, 0].set_xticklabels([])
        axs[0, 1].set_ylabel("rho")
        axs[0, 1].set_ylim([-1.5e-4, 2.5e-4])
        axs[1, 0].plot(x_points, electron_pressure, label="Electrons")
        axs[1, 0].plot(x_points, proton_pressure, label="Protons")
        axs[1, 0].legend(loc='upper right')
        axs[1, 0].set_title(f"Pressure")
        axs[1, 0].set_xlabel("x")
        axs[1, 0].set_ylabel("P")
        axs[1, 0].set_ylim(0, 2.5e-4)
        axs[1, 1].plot(x_points, e_field, label="E")
        axs[1, 1].set_title(f"Electric Field")
        axs[1, 1].set_xlabel("x")
        axs[1, 1].set_ylabel("E")
        axs[1, 1].set_ylim(-1400, 1400)
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/ExpandingGaussian{i:03}.png")
        plt.close(fig)

    txt_data = np.genfromtxt(f"{time_integrator.replace(" ", "_")}/output_lf_0.csv", names=True)
    time = txt_data['Time']
    field_energy = txt_data['Field_Energy']
    fluid_energy = txt_data['Total_Fluid_Energy']
    total_energy = field_energy + fluid_energy
    fluid_kinetic_energy = txt_data['Total_Fluid_Kinetic_Energy']
    fluid_internal_energy = fluid_energy - fluid_kinetic_energy

    fig, ax = plt.subplots()
    ax.plot(time, total_energy, label="Total Energy")
    ax.plot(time, field_energy, label="Field Energy")
    ax.plot(time, fluid_energy, label="Fluid Energy")
    ax.legend()
    ax.set_title(f"Energies Over Time")
    ax.set_xlabel("time")
    ax.set_ylabel("Energy")
    fig.savefig(f"{figures_directory}/Energies.png")
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.plot(time, fluid_energy, label="Total Fluid Energy")
    ax.plot(time, fluid_kinetic_energy, label="Fluid Kinetic Energy")
    ax.plot(time, fluid_internal_energy, label="Fluid Internal Energy")
    ax.legend()
    ax.set_title(f"Fluid Energies Over Time")
    ax.set_xlabel("time")
    ax.set_ylabel("Energy")
    fig.savefig(f"{figures_directory}/FluidEnergies.png")
    plt.close(fig)

if __name__ == "__main__":
    import sys

    time_integrators = ["Forward Euler", "Verlet"]
    # time_integrators = ["Verlet"]
    if "run" in sys.argv[1:]:
        for time_integrator in time_integrators:
            run(sys.argv[2], time_integrator)
    elif "plot" in sys.argv[1:]:
        for time_integrator in time_integrators:
            plot(time_integrator)
    else:
        analyze()
