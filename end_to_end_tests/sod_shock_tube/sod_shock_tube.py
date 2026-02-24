import euler
import euler_exact_riemann_solver
import read_mesh_data
import species

import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.constants import electron_mass
import subprocess

electron_species = species.Species(mass=electron_mass)

mass_density_left = 1.0
pressure_left = 1.0
bulk_velocity_left = np.array([0.0, 0.0, 0.0])

mass_density_right = 0.125
pressure_right = 0.1
bulk_velocity_right = np.array([0.0, 0.0, 0.0])

number_density_left = mass_density_left / electron_species.mass
temperature_left = euler.temperature(number_density_left, pressure_left)

number_density_right = mass_density_right / electron_species.mass
temperature_right = euler.temperature(number_density_right, pressure_right)

sound_speed_left = euler.speed_of_sound(
    electron_species, mass_density_left, pressure_left
)
max_wavespeed_left = np.max(np.abs(bulk_velocity_left)) + sound_speed_left

sound_speed_right = euler.speed_of_sound(
    electron_species, mass_density_right, pressure_right
)
max_wavespeed_right = np.max(np.abs(bulk_velocity_right)) + sound_speed_right

max_wavespeed = max(max_wavespeed_left, max_wavespeed_right)

domain_length = 1.0
num_elements = 100
dx = domain_length / num_elements
cfl = 0.5
dt = cfl * dx / max_wavespeed
initial_final_time = 0.3 * domain_length / max_wavespeed
num_time_steps = int(initial_final_time / dt)
final_time = num_time_steps * dt

discontinuity_location = 0.5 * domain_length

input_deck_contents = f"""
Mesh:
  Type: line
  Lengths: [{domain_length}]
  Number of Elements: [{num_elements}]

Time Stepping:
  Number of Time Steps: {num_time_steps}
  Time Step Size: {dt}

Species:
  neutral_electron:
    Mass: {electron_mass}
    Charge: 0

Euler Fluids:
  Basis Order: 0
  Initial Conditions:
    - Species: [neutral_electron]
      Sod:
        Discontinuity Location: {discontinuity_location}
        Left State:
          Number Density: {number_density_left}
          Temperature: {temperature_left}
          Bulk Velocity: [{bulk_velocity_left[0]}, {bulk_velocity_left[1]}, {bulk_velocity_left[2]}]
        Right State:
          Number Density: {number_density_right}
          Temperature: {temperature_right}
          Bulk Velocity: [{bulk_velocity_right[0]}, {bulk_velocity_right[1]}, {bulk_velocity_right[2]}]
  Boundary Conditions:
    - Side: left
      Type: Reflecting
    - Side: right
      Type: Reflecting

Output:
  Stride: 1

"""


def run(mfpic_executable):
    yaml = "sod_shock_tube.yaml"
    with open(yaml, "w") as input_deck:
        input_deck.write(input_deck_contents)

    result = subprocess.run([mfpic_executable, "-i", yaml])
    result.check_returncode()


def analyze():
    timesteps, mesh_data = read_mesh_data.read_mesh_data()

    points = mesh_data[0]["points"]
    x_points = points[:, 0]

    left_state_primitive = euler.construct_primitive_state(
        number_density_left, bulk_velocity_left, temperature_left
    )
    right_state_primitive = euler.construct_primitive_state(
        number_density_right, bulk_velocity_right, temperature_right
    )
    exact_mass_density, exact_velocity, exact_pressure = (
        euler_exact_riemann_solver.form_exact_solutions(
            left_state_primitive,
            right_state_primitive,
            electron_species,
            discontinuity_location,
        )
    )

    os.makedirs("Figures", exist_ok=True)
    for i in range(len(timesteps)):
        fluid_data = np.transpose(mesh_data[i]["species_0"])
        mass_density = fluid_data[0]
        time = timesteps[i]

        fig, axes = plt.subplots()
        axes.plot(x_points, mass_density, label="Numerical Solution")
        axes.plot(x_points, exact_mass_density(x_points, time), label="Exact Solution")
        axes.legend()
        axes.set_title(f"Mass Density At Time = {time}")
        axes.set_xlabel("x")
        axes.set_ylabel("rho")
        fig.savefig(f"Figures/MassDensity{i:03}.png")
        plt.close(fig)

        bulk_velocity = euler.get_bulk_velocity_from_conservative_state(fluid_data)
        fig, axes = plt.subplots()
        axes.plot(x_points, bulk_velocity[0], label="Numerical Solution")
        axes.plot(x_points, exact_velocity(x_points, time), label="Exact Solution")
        axes.legend()
        axes.set_title(f"Velocity At Time = {time}")
        axes.set_xlabel("x")
        axes.set_ylabel("v")
        fig.savefig(f"Figures/Velocity{i:03}.png")
        plt.close(fig)

        pressure = euler.get_pressure_from_conservative_state(
            fluid_data, electron_species
        )
        fig, axes = plt.subplots()
        axes.plot(x_points, pressure, label="Numerical Solution")
        axes.plot(x_points, exact_pressure(x_points, time), label="Exact Solution")
        axes.set_title(f"Pressure At Time = {time}")
        axes.set_xlabel("x")
        axes.set_ylabel("P")
        fig.savefig(f"Figures/Pressure{i:03}.png")
        plt.close(fig)


if __name__ == "__main__":
    import sys

    if "run" in sys.argv[1:]:
        run(sys.argv[2])
    else:
        analyze()
