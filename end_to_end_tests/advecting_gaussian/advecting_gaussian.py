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

base_num_elements = 25

N2_species = species.Species(mass=4.65e-26, specific_heat_ratio=1.4)
number_density_offset = 1e26
perturbation = 0.001
number_density_height = perturbation * number_density_offset

max_mass_density = (number_density_height + number_density_offset) * N2_species.mass
pressure = 27613
sound_speed = euler.speed_of_sound(N2_species, max_mass_density, pressure)

mach_number = 2.5
velocity = mach_number * sound_speed

max_wavespeed = velocity + sound_speed

final_time = 10 * gaussian_standard_deviation / velocity

basis_order = 0

def exact_mass_density(x, t):
    sigma = gaussian_standard_deviation
    c = gaussian_center
    mask = x < np.fmod(velocity * t, domain_length)
    number_of_periods = int(velocity * t / domain_length)
    shift = (mask + number_of_periods) * domain_length
    exponential = np.exp(-0.5 * np.power(shift + x - velocity * t - c, 2) / np.power(sigma, 2))
    number_density = number_density_height * exponential + number_density_offset
    return N2_species.mass * number_density


def format_mesh_folder_name(refinement_level, time_integrator):
    return f"{time_integrator.replace(" ", "_")}/MeshOutput{refinement_level:02}"


def get_input_deck(refinement_level, time_integrator):
    basis_order = 1 if time_integrator == "SSPERK32" else 0
    max_cfl = 0.3 if time_integrator == "SSPERK32" else 0.95

    num_elements = base_num_elements * refinement_level
    dx = domain_length / num_elements
    dt, num_time_steps = utils.compute_timestepping_that_satisfies_cfl(max_cfl, dx, max_wavespeed, final_time)

    mesh_folder_name = format_mesh_folder_name(refinement_level, time_integrator)

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
  neutral_electron:
    Mass: {N2_species.mass}
    Charge: 0
    Specific Heat Ratio: {N2_species.specific_heat_ratio}

Euler Fluids:
  Basis Order: {basis_order}
  Initial Conditions:
    - Species: [neutral_electron]
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
  Stride: {num_time_steps}
  Mesh Output Folder: {mesh_folder_name}
"""
    return input_deck_contents


def run(mfpic_executable, time_integrator, refinement_levels):
    output_directory = f"{time_integrator.replace(" ", "_")}"
    os.makedirs(output_directory, exist_ok=True)

    for refinement_level in refinement_levels:
        input_deck_contents = get_input_deck(refinement_level, time_integrator)
        yaml = "advecting_gaussian.yaml"
        with open(yaml, "w") as input_deck:
            input_deck.write(input_deck_contents)

        result = subprocess.run([mfpic_executable, "-i", yaml], capture_output=True, text=True)

        log_filename = f"{output_directory}/execute_{refinement_level:02}.log"
        with open(log_filename, "w") as execute_log:
            execute_log.write(result.stdout)
            execute_log.write(result.stderr)

        verification.check_fluid_energy_positive_and_constant('Total_Fluid_Energy')
        verification.check_fluid_energy_positive_and_constant('Total_Fluid_Kinetic_Energy')

        os.rename(yaml, f"{output_directory}/{yaml}_{refinement_level:02}")
        os.rename('output_lf_0.csv', f"{output_directory}/output_lf_0_{refinement_level:02}.csv")
        os.rename(yaml, f"{output_directory}/advecting_gaussian_{refinement_level:02}.yaml")

        result.check_returncode()


def compute_error(data, points, exact_solution):
    numerical_solution = verification.create_1D_interpolater(data, points)
    exact_solution_at_final_time = lambda x: exact_solution(x, final_time)
    error = verification.compute_L1_relative_error_1D(numerical_solution, exact_solution_at_final_time, points)
    return error


def analyze(time_integrator, refinement_levels):
    errors = []
    h_list = []
    for refinement_level in refinement_levels:
        mesh_folder_name = format_mesh_folder_name(refinement_level, time_integrator)
        _, mesh_data = read_mesh_data.read_mesh_data(mesh_folder_name)
        points = mesh_data[-1]["points"]
        x_points = points[:,0]
        dx = x_points[1] - x_points[0]
        h_list.append(dx)

        fluid_data = np.transpose(mesh_data[-1]["species_0_lf_0"])
        error = compute_error(fluid_data[0], x_points, exact_mass_density)
        errors.append(error)

    output_directory = f"{time_integrator.replace(" ", "_")}"
    os.makedirs(output_directory, exist_ok=True)

    rates = verification.compute_convergence_rates(errors, h_list)
    print(f"rates = {rates}")

    rate_data = np.zeros((len(refinement_levels), 4))
    rate_data[:, 0] = refinement_levels
    rate_data[:, 1] = h_list
    rate_data[:, 2] = errors
    rate_data[1:, 3] = rates
    np.savetxt(f"{output_directory}/convergence_rates.txt", rate_data, header = "Refinement_Level H Error Rate")

    expected_convergence_rate = 2 if time_integrator == "SSPERK32" else 1
    figure_name = f"{output_directory}/error_convergence.png"
    verification.plot_errors_and_expected_convergence_rate(h_list, errors, expected_convergence_rate, figure_name)

    tolerance = 0.1
    assert(np.all(rates > expected_convergence_rate - tolerance))


def plot(time_integrator, refinement_levels):
    for refinement_level in refinement_levels:
        mesh_folder_name = format_mesh_folder_name(refinement_level, time_integrator)
        timesteps, mesh_data = read_mesh_data.read_mesh_data(mesh_folder_name)

        points = mesh_data[0]["points"]
        x_points = points[:, 0]
        num_cells = int(0.5 * x_points.shape[0])
        x_plot = np.linspace(0, domain_length, 10 * num_cells)

        figures_directory = f"{time_integrator.replace(" ", "_")}/Figures{refinement_level:02}"
        os.makedirs(figures_directory, exist_ok=True)
        for i, time in enumerate(timesteps):
            fluid_data = np.transpose(mesh_data[i]["species_0_lf_0"])
            mass_density_data = fluid_data[0]

            fig, axes = plt.subplots()
            axes.plot(x_points, mass_density_data, label="Numerical Solution")
            axes.plot(x_plot, exact_mass_density(x_plot, time), label="Exact Solution")
            axes.legend()
            axes.set_title(f"Mass Density At Time = {time}")
            axes.set_xlabel("x")
            axes.set_ylabel("rho")
            fig.savefig(f"{figures_directory}/MassDensity{i:02}.png")
            plt.close(fig)

if __name__ == "__main__":
    import sys

    time_integrators = ["Forward Euler", "Verlet", "SSPERK32"]
    refinement_levels = {}
    refinement_levels["Forward Euler"] = [8, 16, 32]
    refinement_levels["Verlet"] = [8, 16, 32]
    refinement_levels["SSPERK32"] = [2, 4, 8]
    if "run" in sys.argv[1:]:
        mfpic_executable = sys.argv[2]
        for time_integrator in time_integrators:
            run(mfpic_executable, time_integrator, refinement_levels[time_integrator])
    elif "plot" in sys.argv[1:]:
        for time_integrator in time_integrators:
            plot(time_integrator, refinement_levels[time_integrator])
    else:
        for time_integrator in time_integrators:
            analyze(time_integrator, refinement_levels[time_integrator])