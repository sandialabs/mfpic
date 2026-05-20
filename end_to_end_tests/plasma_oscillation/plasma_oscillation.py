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
from scipy.constants import Boltzmann, electron_mass, proton_mass, elementary_charge, epsilon_0
import subprocess

domain_length = 0.01
num_elements = 10
dx = domain_length / num_elements

proton_species = species.Species(charge=elementary_charge, mass=proton_mass)
electron_species = species.Species(charge=-elementary_charge, mass=electron_mass)

number_density = 1e16
temperature = 20
electron_bulk_velocity = 1000
proton_bulk_velocity = -800

electron_mass_density = electron_species.mass * number_density
proton_mass_density = proton_species.mass * number_density
pressure = euler.pressure_from_number_density(number_density, temperature)
electron_sound_speed = euler.speed_of_sound(electron_species, electron_mass_density, pressure)
proton_sound_speed = euler.speed_of_sound(proton_species, proton_mass_density, pressure)

# convenience variables
m0 = electron_species.mass
m1 = proton_species.mass
q0 = electron_species.charge
q1 = proton_species.charge
rho0 = electron_mass_density
rho1 = proton_mass_density
p00 = electron_bulk_velocity * rho0
p10 = proton_bulk_velocity * rho1
eps = epsilon_0


def compute_complex_plasma_frequency():
  denominator = m0 * m1 * eps
  numerator = np.emath.sqrt(-eps * (m1 * m1 * q0 * q0 * rho0 + m0 * m0 * q1 * q1 * rho1))
  complex_plasma_frequency = numerator / denominator
  return complex_plasma_frequency

complex_plasma_frequency = compute_complex_plasma_frequency()
plasma_frequency = np.abs(complex_plasma_frequency)
plasma_period = 2.0 * np.pi / plasma_frequency

final_time = plasma_period
num_time_steps = 20
dt = plasma_period / num_time_steps

electron_cfl = dt * electron_sound_speed / dx
assert(electron_cfl < 1.)

proton_cfl = dt * proton_sound_speed / dx
assert(proton_cfl < 1.)

basis_order = 0


def compute_exact_e_field(t):
  scale = - (m1 * q0 * p00 + m0 * q1 * p10)
  numerator = scale * np.sinh(complex_plasma_frequency * t)
  denominator = np.emath.sqrt(-epsilon_0 * (m1 * m1 * q0 * q0 * rho0 + m0 * m0 * q1 * q1 * rho1))
  e_field = np.real(numerator / denominator)
  return e_field


def compute_exact_p0(t):
  offset = m0 * q1 * (-m1 * p10 * q0 * rho0 + m0 * p00 * q1 * rho1)
  scale = m1 * q0 * (m1 * p00 * q0 + m0 * p10 * q1) * rho0
  numerator = offset + scale * np.cosh(complex_plasma_frequency * t)
  denominator = m1 * m1 * q0 * q0 * rho0 + m0 * m0 * q1 * q1 * rho1
  p0 = np.real(numerator / denominator)
  return p0


def compute_exact_p1(t):
  offset = m1 * q0 * (m1 * p10 * q0 * rho0 - m0 * p00 * q1 * rho1)
  scale = m0 * q1 * (m1 * p00 * q0+ m0 * p10 * q1) * rho1
  numerator = offset + scale * np.cosh(complex_plasma_frequency * t)
  denominator = m1 * m1 * q0 * q0 * rho0 + m0 * m0 * q1 * q1 * rho1
  p1 = np.real(numerator / denominator)
  return p1


input_deck_contents = f"""
Mesh:
  Type: line
  Lengths: [{domain_length}]
  Number of Elements: [{num_elements}]
  Periodic Dimensions: [x]

Time Stepping:
  Number of Time Steps: {num_time_steps}
  Time Step Size: {dt}
  Type: Crank Nicolson

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
    - Species: [electron]
      Constant:
        Number Density: {number_density}
        Temperature: {temperature}
        Bulk Velocity: [{electron_bulk_velocity}, 0., 0.]
    - Species: [proton]
      Constant:
        Number Density: {number_density}
        Temperature: {temperature}
        Bulk Velocity: [{proton_bulk_velocity}, 0., 0.]

Output:
  Stride: 1
"""


def run(mfpic_executable):
    yaml = "plasma_oscillation.yaml"
    with open(yaml, "w") as input_deck:
        input_deck.write(input_deck_contents)

    result = subprocess.run([mfpic_executable, "-i", yaml], capture_output=True, text=True)

    log_filename = f"execute.log"
    with open(log_filename, "w") as execute_log:
        execute_log.write(result.stdout)
        execute_log.write(result.stderr)

    txt_data = np.genfromtxt('output_lf_0.csv', names=True)
    time = txt_data['Time']
    field_energy = txt_data['Field_Energy']
    fluid_energy = txt_data['Total_Fluid_Energy']
    total_energy = field_energy + fluid_energy
    fluid_kinetic_energy = txt_data['Total_Fluid_Kinetic_Energy']
    fluid_internal_energy = fluid_energy - fluid_kinetic_energy

    verification.check_data_positive_and_constant(total_energy)
    verification.check_data_positive_and_constant(fluid_internal_energy)

    result.check_returncode()

def compute_L2_norm_in_time(data):
    L2_norm = np.sqrt(0.5 * np.power(data[0], 2) + np.sum(np.power(data[1:-1], 2)) + 0.5 * np.power(data[-1], 2))
    return L2_norm


def compute_relative_error_in_time(numerical_solution, exact_solution):
    errors = numerical_solution - exact_solution
    absolute_error = compute_L2_norm_in_time(errors)
    exact_norm = compute_L2_norm_in_time(exact_solution)
    relative_error = absolute_error / exact_norm
    return relative_error


def analyze():
    # check error
    mesh_folder_name = "MeshOutput"
    timesteps, mesh_data = read_mesh_data.read_mesh_data(mesh_folder_name)
    e_field = np.array([np.mean(m["e_field_dg_lf_0"][:,0]) for m in mesh_data])
    p0 = np.array([np.mean(m["species_0_lf_0"][:,1]) for m in mesh_data])
    p1 = np.array([np.mean(m["species_1_lf_0"][:,1]) for m in mesh_data])

    exact_e_field = compute_exact_e_field(timesteps)
    exact_p0 = compute_exact_p0(timesteps)
    exact_p1 = compute_exact_p1(timesteps)

    e_field_error = compute_relative_error_in_time(e_field, exact_e_field)
    p0_error = compute_relative_error_in_time(p0, exact_p0)
    p1_error = compute_relative_error_in_time(p1, exact_p1)
    print(f"e_field_error = {e_field_error}")
    print(f"p0_error = {p0_error}")
    print(f"p1_error = {p1_error}")


def plot():
    mesh_folder_name = "MeshOutput"
    timesteps, mesh_data = read_mesh_data.read_mesh_data(mesh_folder_name)

    exact_times = np.linspace(0, final_time)
    exact_e_field = compute_exact_e_field(exact_times)
    exact_p0 = compute_exact_p0(exact_times)
    exact_p1 = compute_exact_p1(exact_times)
    e_field = np.array([np.mean(m["e_field_dg_lf_0"][:,0]) for m in mesh_data])
    p0 = np.array([np.mean(m["species_0_lf_0"][:,1]) for m in mesh_data])
    p1 = np.array([np.mean(m["species_1_lf_0"][:,1]) for m in mesh_data])

    figures_directory = f"Figures"
    os.makedirs(figures_directory, exist_ok=True)

    fig, axes = plt.subplots()
    axes.plot(exact_times, exact_e_field, label="Exact Solution")
    axes.plot(timesteps, e_field, label="Numerical Solution")
    axes.legend()
    axes.set_title(f"Electric Field")
    axes.set_xlabel("t")
    axes.set_ylabel("E")
    fig.savefig(f"{figures_directory}/ElectricField.png")
    plt.close(fig)

    fig, axes = plt.subplots()
    axes.plot(exact_times, exact_p0, label="Exact Solution")
    axes.plot(timesteps, p0, label="Numerical Solution")
    axes.legend()
    axes.set_title(f"Electron Momentum Density")
    axes.set_xlabel("t")
    axes.set_ylabel("p0")
    fig.savefig(f"{figures_directory}/ElectronMomentumDensity.png")
    plt.close(fig)

    fig, axes = plt.subplots()
    axes.plot(exact_times, exact_p1, label="Exact Solution")
    axes.plot(timesteps, p1, label="Numerical Solution")
    axes.legend()
    axes.set_title(f"Proton Momentum Density")
    axes.set_xlabel("t")
    axes.set_ylabel("p1")
    fig.savefig(f"{figures_directory}/ProtonMomentumDensity.png")
    plt.close(fig)

if __name__ == "__main__":
    import sys

    if "run" in sys.argv[1:]:
        run(sys.argv[2])
    if "plot" in sys.argv[1:]:
        plot()
    else:
        analyze()
