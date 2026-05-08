import sys

sys.path.append("../python")

import euler
import euler_exact_riemann_solver
import read_mesh_data
import species
import utils
import verification

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

left_state_primitive = euler.construct_primitive_state(number_density_left, bulk_velocity_left, temperature_left)
right_state_primitive = euler.construct_primitive_state(number_density_right, bulk_velocity_right, temperature_right)

max_wavespeed = euler_exact_riemann_solver.compute_max_wavespeed(left_state_primitive, right_state_primitive, electron_species)

domain_length = 1.0
discontinuity_location = 0.5 * domain_length

max_cfl = 0.8
final_time = 0.4 * domain_length / max_wavespeed

base_num_elements = 50
refinement_levels = [2, 4]


def format_mesh_folder_name(refinement_level, time_integrator):
    return f"{time_integrator.replace(" ", "_")}/MeshOutput{refinement_level:02}"


def get_input_deck(refinement_level, time_integrator):
    num_elements = base_num_elements * refinement_level
    dx = domain_length / num_elements
    dt, num_time_steps = utils.compute_timestepping_that_satisfies_cfl(max_cfl, dx, max_wavespeed, final_time)

    mesh_folder_name = format_mesh_folder_name(refinement_level, time_integrator)

    input_deck_contents = f"""
Mesh:
  Type: line
  Lengths: [{domain_length}]
  Number of Elements: [{num_elements}]

Time Stepping:
  Number of Time Steps: {num_time_steps}
  Time Step Size: {dt}
  Type: {time_integrator}

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
  Mesh Output Folder: {mesh_folder_name}

"""
    return input_deck_contents


def run(mfpic_executable, time_integrator):
    output_directory = f"{time_integrator.replace(" ", "_")}"
    os.makedirs(output_directory, exist_ok=True)

    for refinement_level in refinement_levels:
        input_deck_contents = get_input_deck(refinement_level, time_integrator)
        yaml = "sod_shock_tube.yaml"
        with open(yaml, "w") as input_deck:
            input_deck.write(input_deck_contents)

        result = subprocess.run([mfpic_executable, "-i", yaml], capture_output=True, text=True)

        log_filename = f"{output_directory}/execute_{refinement_level:02}.log"
        with open(log_filename, "w") as execute_log:
            execute_log.write(result.stdout)
            execute_log.write(result.stderr)

        verification.check_fluid_energy_positive_and_constant('Total_Fluid_Energy')
        os.rename('output_lf_0.csv', f"{output_directory}/output_lf_0_{refinement_level:02}.csv")
        os.rename(yaml, f"{output_directory}/sod_shock_tube_{refinement_level:02}.yaml")

        result.check_returncode()


def compute_error(data, points, exact_solution):
    numerical_solution = verification.create_1D_interpolater(data, points)
    exact_solution_at_final_time = lambda x: exact_solution(x, final_time)
    error = verification.compute_L1_relative_error_1D(numerical_solution, exact_solution_at_final_time, points)
    return error


def analyze(time_integrator):
    exact_mass_density, exact_velocity, exact_pressure = euler_exact_riemann_solver.form_exact_solutions(
        left_state_primitive,
        right_state_primitive,
        electron_species,
        discontinuity_location,
    )

    mass_density_errors = []
    velocity_errors = []
    pressure_errors = []
    h_list = []
    for refinement_level in refinement_levels:
        mesh_folder_name = format_mesh_folder_name(refinement_level, time_integrator)
        _, mesh_data = read_mesh_data.read_mesh_data(mesh_folder_name)

        points = mesh_data[-1]["points"]
        x_points = points[:, 0]
        dx = x_points[1] - x_points[0]
        h_list.append(dx)

        fluid_data = np.transpose(mesh_data[-1]["species_0_lf_0"])
        mass_density_error = compute_error(fluid_data[0], x_points, exact_mass_density)

        bulk_velocity_data = euler.get_bulk_velocity_from_conservative_state(fluid_data)
        velocity_error = compute_error(bulk_velocity_data[0], x_points, exact_velocity)

        pressure_data = euler.get_pressure_from_conservative_state(fluid_data, electron_species)
        pressure_error = compute_error(pressure_data, x_points, exact_pressure)

        mass_density_errors.append(mass_density_error)
        velocity_errors.append(velocity_error)
        pressure_errors.append(pressure_error)

    mass_density_rates = verification.compute_convergence_rates(mass_density_errors, h_list)
    velocity_rates = verification.compute_convergence_rates(velocity_errors, h_list)
    pressure_rates = verification.compute_convergence_rates(pressure_errors, h_list)

    print(f"mass_density_rates = {mass_density_rates}")
    print(f"velocity_rates = {velocity_rates}")
    print(f"pressure_rates = {pressure_rates}")

    # since there is a contact wave in this problem the expected convergence rate is
    # only 1/2
    minimum_convergence_rate = 0.5
    assert (
        np.min(mass_density_rates) > minimum_convergence_rate
    ), f"The mass density is not converging at the correct rate. The expected rate is {minimum_convergence_rate} and the actual rates are {mass_density_rates}"
    assert (
        np.min(velocity_rates) > minimum_convergence_rate
    ), f"The velocity is not converging at the correct rate. The expected rate is {minimum_convergence_rate} and the actual rates are {velocity_rates}"
    assert (
        np.min(pressure_rates) > minimum_convergence_rate
    ), f"The pressure is not converging at the correct rate. The expected rate is {minimum_convergence_rate} and the actual rates are {pressure_rates}"


def plot_quantity(data, points, exact_solution, plot_points, time, name, i, figures_directory):
    fig, axes = plt.subplots()
    axes.plot(points, data, label="Numerical Solution")
    axes.plot(plot_points, exact_solution(plot_points, time), label="Exact Solution")
    axes.legend()
    axes.set_title(f"{name} At Time = {time}")
    axes.set_xlabel("x")
    axes.set_ylabel(f"{name}")
    fig.savefig(f"{figures_directory}/{name}{i:03}.png")
    plt.close(fig)


def plot(time_integrator):
    for refinement_level in refinement_levels:
        mesh_folder_name = format_mesh_folder_name(refinement_level, time_integrator)
        timesteps, mesh_data = read_mesh_data.read_mesh_data(mesh_folder_name)

        points = mesh_data[0]["points"]
        x_points = points[:, 0]
        num_cells = int(0.5 * x_points.shape[0])

        exact_mass_density, exact_velocity, exact_pressure = euler_exact_riemann_solver.form_exact_solutions(
            left_state_primitive,
            right_state_primitive,
            electron_species,
            discontinuity_location,
        )
        x_plot = np.linspace(0, domain_length, 10 * num_cells)

        figures_directory = f"{time_integrator.replace(" ", "_")}/Figures{refinement_level:02}"
        os.makedirs(figures_directory, exist_ok=True)
        for i in range(len(timesteps)):
            fluid_data = np.transpose(mesh_data[i]["species_0_lf_0"])
            time = timesteps[i]

            mass_density_data = fluid_data[0]
            plot_quantity(
                mass_density_data,
                x_points,
                exact_mass_density,
                x_plot,
                time,
                "MassDensity",
                i,
                figures_directory,
            )

            bulk_velocity_data = euler.get_bulk_velocity_from_conservative_state(fluid_data)
            plot_quantity(
                bulk_velocity_data[0],
                x_points,
                exact_velocity,
                x_plot,
                time,
                "Velocity",
                i,
                figures_directory,
            )

            pressure_data = euler.get_pressure_from_conservative_state(fluid_data, electron_species)
            plot_quantity(
                pressure_data,
                x_points,
                exact_pressure,
                x_plot,
                time,
                "Pressure",
                i,
                figures_directory,
            )

if __name__ == "__main__":
    import sys

    time_integrators = ["Forward Euler", "Verlet", "SSPERK32"]
    # time_integrators = ["SSPERK32"]
    if "run" in sys.argv[1:]:
        mfpic_executable = sys.argv[2]
        for time_integrator in time_integrators:
            run(mfpic_executable, time_integrator)
    elif "plot" in sys.argv[1:]:
        for time_integrator in time_integrators:
            plot(time_integrator)
    else:
        for time_integrator in time_integrators:
            analyze(time_integrator)