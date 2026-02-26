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

left_state_primitive = euler.construct_primitive_state(
    number_density_left, bulk_velocity_left, temperature_left
)
right_state_primitive = euler.construct_primitive_state(
    number_density_right, bulk_velocity_right, temperature_right
)

max_wavespeed = euler_exact_riemann_solver.compute_max_wavespeed(
    left_state_primitive, right_state_primitive, electron_species
)

domain_length = 1.0
discontinuity_location = 0.5 * domain_length

max_cfl = 0.9
final_time = 0.4 * domain_length / max_wavespeed

def get_input_deck(refinement_level):
    base_num_elements = 50
    num_elements = base_num_elements * refinement_level
    dx = domain_length / num_elements
    max_dt = max_cfl * dx / max_wavespeed

    num_time_steps = int(np.ceil(final_time / max_dt))
    dt = final_time / num_time_steps

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
  Mesh Output Folder: MeshOutput{refinement_level}

"""
    return input_deck_contents


def run(mfpic_executable):

    input_deck_contents = get_input_deck(1)
    yaml = "sod_shock_tube.yaml"
    with open(yaml, "w") as input_deck:
        input_deck.write(input_deck_contents)

    result = subprocess.run([mfpic_executable, "-i", yaml])
    result.check_returncode()


def analyze():
    import verification

    timesteps, mesh_data = read_mesh_data.read_mesh_data()

    final_time = timesteps[-1]

    points = mesh_data[-1]["points"]
    x_points = points[:, 0]

    fluid_data = np.transpose(mesh_data[-1]["species_0"])

    numerical_mass_density = verification.create_1D_interpolater(
        fluid_data[0], x_points
    )
    bulk_velocity_data = euler.get_bulk_velocity_from_conservative_state(fluid_data)
    numerical_velocity = verification.create_1D_interpolater(
        bulk_velocity_data[0], x_points
    )
    pressure_data = euler.get_pressure_from_conservative_state(
        fluid_data, electron_species
    )
    numerical_pressure = verification.create_1D_interpolater(pressure_data, x_points)

    exact_mass_density, exact_velocity, exact_pressure = (
        euler_exact_riemann_solver.form_exact_solutions(
            left_state_primitive,
            right_state_primitive,
            electron_species,
            discontinuity_location,
        )
    )

    exact_mass_density_at_final_time = lambda x: exact_mass_density(x, final_time)
    mass_density_L2_error = verification.compute_L2_error_1D(
        numerical_mass_density, exact_mass_density_at_final_time, x_points
    )
    print(f"mass_density_L2_error = {mass_density_L2_error}")

    exact_velocity_at_final_time = lambda x: exact_velocity(x, final_time)
    velocity_L2_error = verification.compute_L2_error_1D(numerical_velocity, exact_velocity_at_final_time, x_points)
    print(f"velocity_L2_error = {velocity_L2_error}")

    exact_pressure_at_final_time = lambda x: exact_pressure(x, final_time)
    pressure_L2_error = verification.compute_L2_error_1D(numerical_pressure, exact_pressure_at_final_time, x_points)


def plot():
    timesteps, mesh_data = read_mesh_data.read_mesh_data()

    points = mesh_data[0]["points"]
    x_points = points[:, 0]
    num_cells = int(0.5 * x_points.shape[0])

    exact_mass_density, exact_velocity, exact_pressure = (
        euler_exact_riemann_solver.form_exact_solutions(
            left_state_primitive,
            right_state_primitive,
            electron_species,
            discontinuity_location,
        )
    )
    x_plot = np.linspace(0, domain_length, 10 * num_cells)

    os.makedirs("Figures", exist_ok=True)
    for i in range(len(timesteps)):
        fluid_data = np.transpose(mesh_data[i]["species_0"])
        time = timesteps[i]

        mass_density_data = fluid_data[0]
        fig, axes = plt.subplots()
        axes.plot(x_points, mass_density_data, label="Numerical Solution")
        axes.plot(x_plot, exact_mass_density(x_plot, time), label="Exact Solution")
        axes.legend()
        axes.set_title(f"Mass Density At Time = {time}")
        axes.set_xlabel("x")
        axes.set_ylabel("rho")
        fig.savefig(f"Figures/MassDensity{i:03}.png")
        plt.close(fig)

        bulk_velocity_data = euler.get_bulk_velocity_from_conservative_state(fluid_data)
        fig, axes = plt.subplots()
        axes.plot(x_points, bulk_velocity_data[0], label="Numerical Solution")
        axes.plot(x_plot, exact_velocity(x_plot, time), label="Exact Solution")
        axes.legend()
        axes.set_title(f"Velocity At Time = {time}")
        axes.set_xlabel("x")
        axes.set_ylabel("v")
        fig.savefig(f"Figures/Velocity{i:03}.png")
        plt.close(fig)

        pressure_data = euler.get_pressure_from_conservative_state(
            fluid_data, electron_species
        )
        fig, axes = plt.subplots()
        axes.plot(x_points, pressure_data, label="Numerical Solution")
        axes.plot(x_plot, exact_pressure(x_plot, time), label="Exact Solution")
        axes.legend()
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
