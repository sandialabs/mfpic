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

num_elements = 300
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

final_time = 4.e-6

max_cfl = 0.25
max_wavespeed = sound_speed
dt, num_time_steps = utils.compute_timestepping_that_satisfies_cfl(max_cfl, dx, max_wavespeed, final_time)

basis_order = 0

plasma_frequency = np.sqrt(number_density_peak * np.power(electron_species.charge, 2) / (electron_mass * epsilon_0))
plasma_period = 2.0 * np.pi / plasma_frequency

print(f"plasma_period = {plasma_period}")
print(f"dt = {dt}")

num_output_frames = 100
stride = int(np.floor(num_time_steps / num_output_frames))

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
  Stride: {stride}
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
    ssperk32_data = np.genfromtxt('SSPERK32/output_lf_0.csv', names=True)
    verlet_data = np.genfromtxt('Verlet/output_lf_0.csv', names=True)
    cranknicolson_data = np.genfromtxt('Crank_Nicolson/output_lf_0.csv', names=True)

    time = ssperk32_data['Time']

    field_energy_ssperk32 = ssperk32_data['Field_Energy']
    fluid_energy_ssperk32 = ssperk32_data['Total_Fluid_Energy']
    total_energy_ssperk32 = field_energy_ssperk32 + fluid_energy_ssperk32
    fluid_kinetic_energy_ssperk32 = ssperk32_data['Total_Fluid_Kinetic_Energy']
    fluid_internal_energy_ssperk32 = fluid_energy_ssperk32 - fluid_kinetic_energy_ssperk32
    cfl_ssperk32 = ssperk32_data['CFL']

    field_energy_verlet = verlet_data['Field_Energy']
    fluid_energy_verlet = verlet_data['Total_Fluid_Energy']
    total_energy_verlet = field_energy_verlet + fluid_energy_verlet
    fluid_kinetic_energy_verlet = verlet_data['Total_Fluid_Kinetic_Energy']
    fluid_internal_energy_verlet = fluid_energy_verlet - fluid_kinetic_energy_verlet
    cfl_verlet = verlet_data['CFL']

    field_energy_cranknicolson = cranknicolson_data['Field_Energy']
    fluid_energy_cranknicolson = cranknicolson_data['Total_Fluid_Energy']
    total_energy_cranknicolson = field_energy_cranknicolson + fluid_energy_cranknicolson
    fluid_kinetic_energy_cranknicolson = cranknicolson_data['Total_Fluid_Kinetic_Energy']
    fluid_internal_energy_cranknicolson = fluid_energy_cranknicolson - fluid_kinetic_energy_cranknicolson
    cfl_cranknicolson = cranknicolson_data['CFL']

    figures_directory = "Figures"
    os.makedirs(figures_directory, exist_ok=True)

    fig, ax = plt.subplots()
    ax.plot(time, total_energy_ssperk32, label='SSPERK32')
    ax.plot(time, total_energy_verlet, label='Verlet')
    ax.plot(time, total_energy_cranknicolson, label='Crank Nicolson')
    ax.legend()
    ax.set_title(f"Total Energy")
    ax.set_xlabel("time")
    ax.set_ylabel("Energy")
    fig.tight_layout()
    fig.savefig(f"{figures_directory}/TotalEnergy.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.plot(time[1:], cfl_ssperk32[1:], label='SSPERK32')
    ax.plot(time[1:], cfl_verlet[1:], label='Verlet')
    ax.plot(time[1:], cfl_cranknicolson[1:], label='Crank Nicolson')
    ax.hlines(1., 0., 1., transform=ax.get_yaxis_transform(), colors='k', linestyles='dashed')
    ax.legend()
    ax.set_title(f"CFL")
    ax.set_xlabel("time")
    ax.set_ylabel("CFL")
    fig.tight_layout()
    fig.savefig(f"{figures_directory}/CFL.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.plot(time, fluid_internal_energy_ssperk32, label="SSPERK32")
    ax.plot(time, fluid_internal_energy_verlet, label="Verlet")
    ax.plot(time, fluid_internal_energy_cranknicolson, label="Crank Nicolson")
    ax.legend()
    ax.set_title(f"Fluid Internal Energy over Time")
    ax.set_xlabel("time")
    ax.set_ylabel("Internal Energy")
    fig.tight_layout()
    fig.savefig(f"{figures_directory}/InternalEnergy.pdf")
    plt.close(fig)


def plot(time_integrator):
    mesh_folder_name = format_mesh_folder_name(time_integrator)
    timesteps, mesh_data = read_mesh_data.read_mesh_data(mesh_folder_name)

    points = mesh_data[0]["points"]
    x_points = points[:, 0]

    all_electron_data = np.array([np.transpose(m["species_0_lf_0"]) for m in mesh_data])
    all_proton_data = np.array([np.transpose(m["species_1_lf_0"]) for m in mesh_data])

    max_electron_number_density = np.max([np.max(euler.get_number_density_from_conservative_state(ed, electron_species)) for ed in all_electron_data])
    max_proton_number_density = np.max([np.max(euler.get_number_density_from_conservative_state(d, proton_species)) for d in all_proton_data])
    max_number_density = max(max_electron_number_density, max_proton_number_density)

    min_electron_bulk_velocity = np.min([np.min(euler.get_bulk_velocity_from_conservative_state(d)) for d in all_electron_data])
    max_electron_bulk_velocity = np.max([np.max(euler.get_bulk_velocity_from_conservative_state(d)) for d in all_electron_data])

    max_electron_pressure = np.max([np.max(euler.get_pressure_from_conservative_state(d, electron_species)) for d in all_electron_data])
    min_electron_pressure = np.min([np.min(euler.get_pressure_from_conservative_state(d, electron_species)) for d in all_electron_data])
    max_proton_pressure = np.max([np.max(euler.get_pressure_from_conservative_state(d, proton_species)) for d in all_proton_data])
    min_proton_pressure = np.min([np.min(euler.get_pressure_from_conservative_state(d, proton_species)) for d in all_proton_data])
    max_pressure = max(max_electron_pressure, max_proton_pressure)
    min_pressure = min(min_electron_pressure, min_proton_pressure, 0)
    print(f"min_pressure = {min_pressure}")

    max_electron_internal_energy_density = np.max([np.max(euler.get_internal_energy_density_from_conservative_state(d)) for d in all_electron_data])
    min_electron_internal_energy_density = np.min([np.min(euler.get_internal_energy_density_from_conservative_state(d)) for d in all_electron_data])
    max_proton_internal_energy_density = np.max([np.max(euler.get_internal_energy_density_from_conservative_state(d)) for d in all_proton_data])
    min_proton_internal_energy_density = np.min([np.min(euler.get_internal_energy_density_from_conservative_state(d)) for d in all_proton_data])
    max_internal_energy_density = max(max_electron_internal_energy_density, max_proton_internal_energy_density)
    min_internal_energy_density = min(min_electron_internal_energy_density, min_proton_internal_energy_density, 0)
    print(f"min_internal_energy_density = {min_internal_energy_density}")

    max_electron_temperature = np.max([np.max(euler.get_temperature_from_conservative_state(d, electron_species)) for d in all_electron_data])
    min_electron_temperature = np.min([np.min(euler.get_temperature_from_conservative_state(d, electron_species)) for d in all_electron_data])
    max_proton_temperature = np.max([np.max(euler.get_temperature_from_conservative_state(d, proton_species)) for d in all_proton_data])
    min_proton_temperature = np.min([np.min(euler.get_temperature_from_conservative_state(d, proton_species)) for d in all_proton_data])
    max_temperature = max(max_electron_temperature, max_proton_temperature)
    min_temperature = min(min_electron_temperature, min_proton_temperature, 0)
    print(f"min_temperature = {min_temperature}")

    max_charge_density = np.max([np.max(euler.get_number_density_from_conservative_state(pd, proton_species) - euler.get_number_density_from_conservative_state(ed, electron_species)) for (ed, pd) in zip(all_electron_data, all_proton_data)]) * elementary_charge
    min_charge_density = np.min([np.min(euler.get_number_density_from_conservative_state(pd, proton_species) - euler.get_number_density_from_conservative_state(ed, electron_species)) for (ed, pd) in zip(all_electron_data, all_proton_data)]) * elementary_charge

    max_potential = np.max([np.max(m["electrostatic_potential_lf_0"]) for m in mesh_data])
    min_potential = np.min([np.min(m["electrostatic_potential_lf_0"]) for m in mesh_data])
    if max_potential == min_potential:
        min_potential = -1
        max_potential = 1

    max_e_field = np.max([np.max(m['E_0_lf_0']) for m in mesh_data])
    min_e_field = np.min([np.min(m['E_0_lf_0']) for m in mesh_data])
    if max_e_field == min_e_field:
        min_e_field = -1
        max_e_field = 1

    max_e_field_dg = np.max([np.max(m['e_field_dg_lf_0'][:,0]) for m in mesh_data])
    min_e_field_dg = np.min([np.min(m['e_field_dg_lf_0'][:,0]) for m in mesh_data])
    if max_e_field_dg == min_e_field_dg:
      min_e_field_dg = -1
      max_e_field_dg = 1

    max_ghost_charge_density = np.max([np.max(m['ghost_charge_density_lf_0']) for m in mesh_data])
    min_ghost_charge_density = np.min([np.min(m['ghost_charge_density_lf_0']) for m in mesh_data])
    max_ghost_charge_density_2 = np.max([np.max(m['ghost_charge_density_2_lf_0']) for m in mesh_data])
    min_ghost_charge_density_2 = np.min([np.min(m['ghost_charge_density_2_lf_0']) for m in mesh_data])

    fontdict = {'weight': 'bold', 'size': 12}

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
        electron_temperature = euler.get_temperature_from_conservative_state(electron_data, electron_species)
        proton_temperature = euler.get_temperature_from_conservative_state(proton_data, proton_species)

        electron_number_density = electron_mass_density / electron_species.mass
        proton_number_density = proton_mass_density / proton_species.mass
        charge_density = (proton_number_density - electron_number_density) * elementary_charge

        potential = mesh_data[i]["electrostatic_potential_lf_0"]
        e_field = mesh_data[i]["E_0_lf_0"]
        e_field_dg = mesh_data[i]["e_field_dg_lf_0"][:,0]

        ghost_charge_density = mesh_data[i]["ghost_charge_density_lf_0"]
        ghost_charge_density_2 = mesh_data[i]["ghost_charge_density_2_lf_0"]

        fig, ax = plt.subplots()
        ax.plot(x_points, electron_number_density, label="Electrons")
        ax.plot(x_points, proton_number_density, label="Protons")
        ax.legend(loc='upper right')
        ax.set_title(f"Number Density")
        ax.set_xlabel("x")
        ax.set_ylabel("n")
        ax.set_ylim(0, 1.1 * max_number_density )
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/NumberDensity{i:03}.pdf")
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.plot(x_points, electron_bulk_velocity[0], label="Electrons")
        ax.plot(x_points, proton_bulk_velocity[0], label="Protons")
        ax.legend(loc='upper right')
        ax.set_title(f"Bulk Velocity")
        ax.set_xlabel("x")
        ax.set_ylabel("v")
        ax.set_ylim(1.1 * min_electron_bulk_velocity, 1.1 * max_electron_bulk_velocity)
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/BulkVelocity{i:03}.pdf")
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.plot(x_points, electron_pressure, label="Electrons")
        ax.plot(x_points, proton_pressure, label="Protons")
        ax.legend(loc='upper right')
        ax.set_title(f"Pressure")
        ax.set_xlabel("x")
        ax.set_ylabel("P")
        ax.set_ylim(1.1 * min_pressure, 1.1 * max_pressure)
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/Pressure{i:03}.pdf")
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.plot(x_points, electron_internal_energy_density, label="Electrons")
        ax.plot(x_points, proton_internal_energy_density, label="Protons")
        ax.legend(loc='upper right')
        ax.set_title(f"Internal Energy Density")
        ax.set_xlabel("x")
        ax.set_ylabel("e")
        ax.set_ylim(1.1 * min_internal_energy_density, 1.1 * max_internal_energy_density)
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/InternalEnergyDensity{i:03}.pdf")
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.plot(x_points, electron_temperature, label="Electrons")
        ax.plot(x_points, proton_temperature, label="Protons")
        ax.legend(loc='upper right')
        ax.set_title(f"Temperature")
        ax.set_xlabel("x")
        ax.set_ylabel("T")
        ax.set_ylim(1.1 * min_temperature, 1.1 * max_temperature)
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/Temperature{i:03}.pdf")
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.plot(x_points, charge_density, label="Charge Density")
        # ax.legend()
        ax.set_title(f"Charge Density")
        ax.set_xlabel("x")
        ax.set_ylabel("rho")
        ax.set_ylim([1.1 * min_charge_density, 1.1 * max_charge_density])
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/ChargeDensity{i:03}.pdf")
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.plot(x_points, potential, label="Potential")
        # ax.legend()
        ax.set_title(f"Potential")
        ax.set_xlabel("x")
        ax.set_ylabel("Phi")
        ax.set_ylim(1.1 * min_potential, 1.1 * max_potential)
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/Potential{i:03}.pdf")
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.plot(x_points, e_field, label="E")
        # ax.legend()
        ax.set_title(f"Electric Field")
        ax.set_xlabel("x")
        ax.set_ylabel("E")
        ax.set_ylim(1.1 * min_e_field, 1.1 * max_e_field)
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/ElectricField{i:03}.pdf")
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.plot(x_points, e_field_dg, label="E")
        # ax.legend()
        ax.set_title(f"Electric Field")
        ax.set_xlabel("x")
        ax.set_ylabel("E")
        ax.set_ylim(1.1 * min_e_field_dg, 1.1 * max_e_field_dg)
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/DGElectricField{i:03}.pdf")
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.plot(x_points, ghost_charge_density, label="Gauss' Law Error")
        # ax.legend()
        ax.set_title(f"Error in Gauss' Law")
        ax.set_xlabel("x")
        ax.set_ylabel("charge density")
        ax.set_ylim(1.1 * min_ghost_charge_density, 1.1 * max_ghost_charge_density)
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/GhostChargeDensityFromPotential{i:03}.pdf")
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.plot(x_points, ghost_charge_density_2, label="Gauss' Law Error")
        # ax.legend()
        ax.set_title(f"Error in Gauss' Law")
        ax.set_xlabel("x")
        ax.set_ylabel("charge density")
        ax.set_ylim(1.1 * min_ghost_charge_density_2, 1.1 * max_ghost_charge_density_2)
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/GhostChargeDensityFromEField{i:03}.pdf")
        plt.close(fig)

        fig, axs = plt.subplots(2, 2, figsize=(16, 9))
        axs[0, 0].plot(x_points, electron_number_density, label="Electrons")
        axs[0, 0].plot(x_points, proton_number_density, label="Protons")
        axs[0, 0].legend(loc='upper right', prop=fontdict)
        axs[0, 0].set_title(f"Number Density", fontdict=fontdict)
        axs[0, 0].set_xticklabels([])
        axs[0, 0].set_ylabel("n", fontdict=fontdict)
        axs[0, 0].set_ylim(0, 1.1 * max_number_density)
        axs[0, 0].tick_params(labelsize='medium')
        axs[0, 1].plot(x_points, charge_density, label="Charge Density")
        axs[0, 1].set_title(f"Charge Density", fontdict=fontdict)
        axs[0, 1].set_xticklabels([])
        axs[0, 1].set_ylabel("rho", fontdict=fontdict)
        axs[0, 1].set_ylim([1.1 * min_charge_density, 1.1 * max_charge_density])
        axs[0, 1].tick_params(labelsize='medium')
        axs[1, 0].plot(x_points, electron_internal_energy_density, label="Electrons")
        axs[1, 0].plot(x_points, proton_internal_energy_density, label="Protons")
        axs[1, 0].legend(loc='upper right', prop=fontdict)
        axs[1, 0].set_title(f"Internal Energy Density", fontdict=fontdict)
        axs[1, 0].set_xlabel("x", fontdict=fontdict)
        axs[1, 0].set_ylabel("e", fontdict=fontdict)
        axs[1, 0].set_ylim(0, 1.1 * max_internal_energy_density)
        axs[1, 0].tick_params(labelsize='medium')
        if time_integrator == "Crank Nicolson":
            axs[1, 1].plot(x_points, e_field_dg, label="E")
            axs[1, 1].set_ylim(1.1 * min_e_field_dg, 1.1 * max_e_field_dg)
        else:
            axs[1, 1].plot(x_points, e_field, label="E")
            axs[1, 1].set_ylim(1.1 * min_e_field, 1.1 * max_e_field)
        axs[1, 1].set_title(f"Electric Field", fontdict=fontdict)
        axs[1, 1].set_xlabel("x", fontdict=fontdict)
        axs[1, 1].set_ylabel("E", fontdict=fontdict)
        axs[1, 1].tick_params(labelsize='medium')
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/ExpandingGaussian{i:03}.pdf")
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
    fig.tight_layout()
    fig.savefig(f"{figures_directory}/Energies.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.plot(time, fluid_energy, label="Total Fluid Energy")
    ax.plot(time, fluid_kinetic_energy, label="Fluid Kinetic Energy")
    ax.plot(time, fluid_internal_energy, label="Fluid Internal Energy")
    ax.legend()
    ax.set_title(f"Fluid Energies Over Time")
    ax.set_xlabel("time")
    ax.set_ylabel("Energy")
    fig.tight_layout()
    fig.savefig(f"{figures_directory}/FluidEnergies.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.plot(time, total_energy, label="Total Energy")
    ax.plot(time, fluid_internal_energy, label="Fluid Internal Energy")
    ax.legend()
    ax.set_title(f"Energies Over Time")
    ax.set_xlabel("time")
    ax.set_ylabel("Energy")
    fig.tight_layout()
    fig.savefig(f"{figures_directory}/TotalEnergyAndInternalEnergy.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.semilogy(time, total_energy)
    ax.set_title("Total Energy")
    ax.set_xlabel("time")
    ax.set_ylabel("Total Energy")
    fig.tight_layout()
    fig.savefig(f"{figures_directory}/TotalEnergy.pdf")
    plt.close(fig)


def plot_presentation(time_integrator):
    txt_data = np.genfromtxt(f"{time_integrator.replace(" ", "_")}/output_lf_0.csv", names=True)

    figures_directory = f"{time_integrator.replace(" ", "_")}/PresentationFigures"
    os.makedirs(figures_directory, exist_ok=True)

    time = txt_data['Time']
    field_energy = txt_data['Field_Energy']
    fluid_energy = txt_data['Total_Fluid_Energy']
    total_energy = field_energy + fluid_energy
    fluid_kinetic_energy = txt_data['Total_Fluid_Kinetic_Energy']
    fluid_internal_energy = fluid_energy - fluid_kinetic_energy
    cfl = txt_data['CFL']

    fontdict = {'weight': 'bold', 'size': 12}

    fig, ax = plt.subplots()
    ax.plot(time, total_energy, label="Total Energy", linewidth=2)
    ax.plot(time, fluid_internal_energy, label="Fluid Internal Energy", linewidth=2)
    ax.legend()
    ax.set_title(f"Energies Over Time", fontdict=fontdict)
    ax.set_xlabel("Time", fontdict=fontdict)
    ax.set_ylabel("Energy", fontdict=fontdict)
    ax.tick_params(labelsize='medium')
    fig.tight_layout()
    fig.savefig(f"{figures_directory}/TotalEnergyAndInternalEnergy.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.semilogy(time, total_energy, label="Total Energy", linewidth=2)
    ax.set_title(f"Total Energy Over Time", fontdict=fontdict)
    ax.set_xlabel("Time", fontdict=fontdict)
    ax.set_ylabel("Energy", fontdict=fontdict)
    ax.tick_params(labelsize='medium')
    ax.set_ylim(0.99 * np.min(total_energy), 1.01 * np.max(total_energy))
    fig.tight_layout()
    fig.savefig(f"{figures_directory}/TotalEnergy.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.plot(time[1:], cfl[1:], label="CFL", linewidth=2)
    ax.hlines(1., 0., 1., transform=ax.get_yaxis_transform(), colors='k', linestyles='dashed')
    # ax.legend()
    ax.set_title(f"CFL Number Over Time", fontdict=fontdict)
    ax.set_xlabel("Time", fontdict=fontdict)
    ax.set_ylabel("CFL", fontdict=fontdict)
    ax.tick_params(labelsize='medium')
    fig.tight_layout()
    fig.savefig(f"{figures_directory}/CFL.pdf")
    plt.close(fig)

    mesh_folder_name = format_mesh_folder_name(time_integrator)
    timesteps, mesh_data = read_mesh_data.read_mesh_data(mesh_folder_name)

    points = mesh_data[0]["points"]
    x_points = points[:, 0]

    electron_data = np.transpose(mesh_data[-1]["species_0_lf_0"])
    proton_data = np.transpose(mesh_data[-1]["species_1_lf_0"])

    electron_internal_energy_density = euler.get_internal_energy_density_from_conservative_state(electron_data)
    proton_internal_energy_density = euler.get_internal_energy_density_from_conservative_state(proton_data)

    fig, ax = plt.subplots()
    ax.plot(x_points, electron_internal_energy_density, label="Electrons", linewidth=2)
    ax.plot(x_points, proton_internal_energy_density, label="Protons", linewidth=2)

    ax.legend(loc='upper right')
    ax.set_title(f"Internal Energy Density", fontdict=fontdict)
    ax.set_xlabel("x", fontdict=fontdict)
    ax.set_ylabel("e", fontdict=fontdict)
    ax.hlines(0., 0., 1., transform=ax.get_yaxis_transform(), colors='k', linestyles='dashed')
    ax.tick_params(labelsize='medium')
    fig.tight_layout()
    fig.savefig(f"{figures_directory}/InternalEnergyDensity.pdf")
    plt.close(fig)

def plot_comparisons():
    txt_data_ssperk32 = np.genfromtxt("SSPERK32/output_lf_0.csv", names=True)
    txt_data_verlet = np.genfromtxt("Verlet/output_lf_0.csv", names=True)

    time = txt_data_ssperk32['Time']
    field_energy_ssperk32 = txt_data_ssperk32['Field_Energy']
    fluid_energy_ssperk32 = txt_data_ssperk32['Total_Fluid_Energy']
    total_energy_ssperk32 = field_energy_ssperk32 + fluid_energy_ssperk32
    fluid_kinetic_energy_ssperk32 = txt_data_ssperk32['Total_Fluid_Kinetic_Energy']
    fluid_internal_energy_ssperk32 = fluid_energy_ssperk32 - fluid_kinetic_energy_ssperk32

    field_energy_verlet = txt_data_verlet['Field_Energy']
    fluid_energy_verlet = txt_data_verlet['Total_Fluid_Energy']
    total_energy_verlet = field_energy_verlet + fluid_energy_verlet
    fluid_kinetic_energy_verlet = txt_data_verlet['Total_Fluid_Kinetic_Energy']
    fluid_internal_energy_verlet = fluid_energy_verlet - fluid_kinetic_energy_verlet

    figures_directory = "Figures"
    os.makedirs(figures_directory, exist_ok=True)

    fontdict = {'weight': 'bold', 'size': 12}

    fig, ax = plt.subplots()
    ax.plot(time, total_energy_ssperk32, label="Total Energy (SSPERK32)", linewidth=2)
    ax.plot(time, fluid_internal_energy_ssperk32, label="Fluid Internal Energy (SSPERK32)", linewidth=2)
    ax.plot(time, total_energy_verlet, label="Total Energy (Verlet)", linewidth=2)
    ax.plot(time, fluid_internal_energy_verlet, label="Fluid Internal Energy (Verlet)", linewidth=2)
    ax.legend()
    ax.set_title(f"Energies Over Time", fontdict=fontdict)
    ax.set_xlabel("Time", fontdict=fontdict)
    ax.set_ylabel("Energy", fontdict=fontdict)
    ax.tick_params(labelsize='medium')
    fig.tight_layout()
    fig.savefig(f"{figures_directory}/TotalEnergyAndInternalEnergySSPERK32AndVerlet.pdf")
    plt.close(fig)

    mesh_folder_name_cn = format_mesh_folder_name("Crank Nicolson")
    timesteps, mesh_data_cn = read_mesh_data.read_mesh_data(mesh_folder_name_cn)
    mesh_folder_name_verlet = format_mesh_folder_name("Verlet")
    _, mesh_data_verlet = read_mesh_data.read_mesh_data(mesh_folder_name_verlet)
    mesh_folder_name_ssperk32 = format_mesh_folder_name("SSPERK32")
    _, mesh_data_ssperk32 = read_mesh_data.read_mesh_data(mesh_folder_name_ssperk32)

    points = mesh_data_cn[0]["points"]
    x_points = points[:, 0]

    i = 20
    ghost_charge_density_cn = mesh_data_cn[i]["ghost_charge_density_2_lf_0"]
    ghost_charge_density_verlet = mesh_data_verlet[i]["ghost_charge_density_lf_0"]

    fig, ax = plt.subplots()
    ax.semilogy(x_points, ghost_charge_density_verlet, label="Verlet")
    ax.semilogy(x_points, ghost_charge_density_cn, label="Crank Nicolson")
    ax.legend(prop=fontdict)
    ax.set_title(f"Error in Gauss' Law", fontdict=fontdict)
    ax.set_xlabel("x", fontdict=fontdict)
    ax.set_ylabel("ghost charge density", fontdict=fontdict)
    ax.set_yscale('symlog', linthresh=1e-18)
    fig.tight_layout()
    fig.savefig(f"{figures_directory}/GhostChargeDensityCrankNicolsonAndVerlet.pdf")
    plt.close(fig)

    electron_data_cn = np.transpose(mesh_data_cn[i]["species_0_lf_0"])
    proton_data_cn = np.transpose(mesh_data_cn[i]["species_1_lf_0"])
    electron_mass_density_cn = electron_data_cn[0]
    proton_mass_density_cn = proton_data_cn[0]
    electron_number_density_cn = electron_mass_density_cn / electron_species.mass
    proton_number_density_cn = proton_mass_density_cn / proton_species.mass
    charge_density_cn = (proton_number_density_cn - electron_number_density_cn) * elementary_charge
    weak_div_E_cn = mesh_data_cn[i]["weak_div_E_lf_0"]

    fig, ax = plt.subplots()
    ax.plot(x_points, ghost_charge_density_cn, label="Ghost Charge Density")
    ax.plot(x_points, charge_density_cn, label="Charge Density")
    ax.plot(x_points, weak_div_E_cn, label="Div(D)")
    ax.legend(prop={'weight': 'bold', 'size': 10})
    ax.set_title(f"Error in Gauss' Law", fontdict=fontdict)
    ax.set_xlabel("x", fontdict=fontdict)
    ax.set_ylabel("charge density", fontdict=fontdict)
    fig.tight_layout()
    fig.savefig(f"{figures_directory}/GhostChargeDensityCrankNicolson.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.plot(x_points, ghost_charge_density_cn, label="Ghost Charge Density")
    ax.plot(x_points, charge_density_cn, label="Charge Density")
    ax.plot(x_points, weak_div_E_cn, label="Div(D)")
    ax.legend(prop={'weight': 'bold', 'size': 10})
    ax.set_title(f"Error in Gauss' Law", fontdict=fontdict)
    ax.set_xlabel("x", fontdict=fontdict)
    ax.set_ylabel("charge density", fontdict=fontdict)
    ax.set_yscale('symlog', linthresh=1e-8)
    fig.tight_layout()
    fig.savefig(f"{figures_directory}/GhostChargeDensityCrankNicolsonLog.pdf")
    plt.close(fig)

    ghost_charge_density_ssperk32 = mesh_data_ssperk32[i]['ghost_charge_density_lf_0']
    fig, ax = plt.subplots()
    ax.plot(x_points, ghost_charge_density_cn, label="Crank Nicolson")
    ax.plot(x_points, ghost_charge_density_ssperk32, label="SSPERK32")
    ax.legend(prop={'weight': 'bold', 'size': 10})
    ax.set_title(f"Error in Gauss' Law", fontdict=fontdict)
    ax.set_xlabel("x", fontdict=fontdict)
    ax.set_ylabel("ghost charge density", fontdict=fontdict)
    ax.set_yscale('symlog', linthresh=1e-18)
    fig.tight_layout()
    fig.savefig(f"{figures_directory}/GhostChargeDensitySSPERK32AndCrankNicolson.pdf")
    plt.close(fig)

    for i in range(21):
        electron_data_verlet = np.transpose(mesh_data_verlet[i]["species_0_lf_0"])
        proton_data_verlet = np.transpose(mesh_data_verlet[i]["species_1_lf_0"])
        electron_mass_density_verlet = electron_data_verlet[0]
        proton_mass_density_verlet = proton_data_verlet[0]
        electron_number_density_verlet = electron_mass_density_verlet / electron_species.mass
        proton_number_density_verlet = proton_mass_density_verlet / proton_species.mass

        electron_data_cn = np.transpose(mesh_data_cn[i]["species_0_lf_0"])
        proton_data_cn = np.transpose(mesh_data_cn[i]["species_1_lf_0"])
        electron_mass_density_cn = electron_data_cn[0]
        proton_mass_density_cn = proton_data_cn[0]
        electron_number_density_cn = electron_mass_density_cn / electron_species.mass
        proton_number_density_cn = proton_mass_density_cn / proton_species.mass

        e_field_cn = mesh_data_cn[i]["e_field_dg_lf_0"][:,0]

        electron_data_ssperk32 = np.transpose(mesh_data_ssperk32[i]["species_0_lf_0"])
        proton_data_ssperk32 = np.transpose(mesh_data_ssperk32[i]["species_1_lf_0"])
        electron_mass_density_ssperk32 = electron_data_ssperk32[0]
        proton_mass_density_ssperk32 = proton_data_ssperk32[0]
        electron_number_density_ssperk32 = electron_mass_density_ssperk32 / electron_species.mass
        proton_number_density_ssperk32 = proton_mass_density_ssperk32 / proton_species.mass

        e_field_ssperk32 = mesh_data_ssperk32[i]["E_0_lf_0"]

        fig, ax = plt.subplots(figsize=(8.8, 6.6))
        ax.plot(x_points, electron_number_density_cn, label="Electrons(Crank Nicolson)")
        ax.plot(x_points, proton_number_density_cn, label="Protons(Crank Nicolson)")
        ax.plot(x_points, electron_number_density_ssperk32, label="Electrons(SSPERK32)")
        ax.plot(x_points, proton_number_density_ssperk32, label="Protons(SSPERK32)")
        ax.legend(loc='upper left', prop={'weight':'bold'})
        ax.set_title(f"Number Density", fontdict=fontdict)
        ax.set_xlabel("x", fontdict=fontdict)
        ax.set_ylabel("n", fontdict=fontdict)
        ax.set_ylim(0, 1e16)
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/NumberDensitySSPERK32AndCrankNicolson{i:02}.pdf")
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.plot(x_points, electron_number_density_cn, label="Electrons(Crank Nicolson)")
        ax.plot(x_points, proton_number_density_cn, label="Protons(Crank Nicolson)")
        ax.plot(x_points, electron_number_density_verlet, label="Electrons(Verlet)")
        ax.plot(x_points, proton_number_density_verlet, label="Protons(Verlet)")
        ax.legend(loc='upper left', prop={'weight':'bold'})
        ax.set_title(f"Number Density", fontdict=fontdict)
        ax.set_xlabel("x", fontdict=fontdict)
        ax.set_ylabel("n", fontdict=fontdict)
        ax.set_ylim(0, 1e16)
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/NumberDensityVerletAndCrankNicolson{i:02}.pdf")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8.8, 6.6))
        ax.plot(x_points, e_field_cn, label="Crank Nicolson")
        ax.plot(x_points, e_field_ssperk32, label="SSPERK32")
        ax.legend(loc='upper left', prop={'weight':'bold'})
        ax.set_title(f"Electric Field", fontdict=fontdict)
        ax.set_xlabel("x", fontdict=fontdict)
        ax.set_ylabel("E", fontdict=fontdict)
        # ax.set_ylim(0, 1e16)
        fig.tight_layout()
        fig.savefig(f"{figures_directory}/EFieldSSPERK32AndCrankNicolson{i:02}.pdf")
        plt.close(fig)

if __name__ == "__main__":
    import sys

    time_integrators = ["Crank Nicolson", "SSPERK32", "Verlet"]
    # time_integrators = ["Verlet"]
    # time_integrators = ["Forward Euler"]
    if "run" in sys.argv[1:]:
        for time_integrator in time_integrators:
            run(sys.argv[2], time_integrator)
    elif "plot" in sys.argv[1:]:
        for time_integrator in time_integrators:
            plot(time_integrator)
    elif "plot_presentation" in sys.argv[1:]:
        plot_comparisons()
        for time_integrator in time_integrators:
            plot_presentation(time_integrator)
    # else:
    #     analyze()
