import sys

sys.path.append("../python")

import euler
import euler_exact_riemann_solver
import read_mesh_data
import species
import utils
import verification

import itertools
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.constants import electron_mass
import subprocess
import vrstudy

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

num_elements = 100

dt_times_collision_frequencies = [0.1]

particles_per_cell = 100

def format_test_name(dt_times_collision_frequency : float, use_variance_reduction : bool):
    return f'nu_x_dt_{dt_times_collision_frequency}_vr_{use_variance_reduction}'


def get_input_deck(dt_times_collision_frequency : float, use_variance_reduction : bool):
    dx = domain_length / num_elements
    dt, num_time_steps = utils.compute_timestepping_that_satisfies_cfl(max_cfl, dx, max_wavespeed, final_time)
    collision_frequency = dt_times_collision_frequency / dt

    no_variance_reduction = f"""
Variance Reduction:
  Strategy: none
    """

    with_variance_reduction = f"""
Variance Reduction:
  Strategy: Euler Fluid
    """

    collision_block = f"""
Collisions:
    - BGK Relaxation:
        Species: neutral_electron
        Collision Frequency: {collision_frequency}
    """

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

Particles:
  Default Boundary Condition: Reflecting
  Initial Conditions:
    - Species: [neutral_electron]
      Number of Macroparticles per Species: {particles_per_cell * num_elements}
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

{collision_block if collision_frequency > 0 else ""}

Output:
  Stride: 1
  Mesh Output Folder: MeshOutput_{format_test_name(dt_times_collision_frequency, use_variance_reduction)}

{with_variance_reduction if use_variance_reduction else no_variance_reduction}

"""
    return input_deck_contents


def run(mfpic_executable):
    output_directory = "./"
    os.makedirs(output_directory, exist_ok=True)

    study_toml_contents = f"""
[study]
name = "neutral_electron"
mfpic = "{mfpic_executable}"              # executable name on PATH, or an absolute path
members = 112                # ensemble size per case
base_seed = 1000             # member i uses Particles: Seed = base_seed + i in every case
jobs = 112                   # concurrent mfpic runs (each pinned to one core, OMP_NUM_THREADS=1)
pin_cores = true
species = ["neutral_electron"]
keep_particle_dumps = [1]    # members that keep particles.h5part (all others turn it off)
keep_mesh_output = true

[collect]
# MeshOutput point fields to keep (the matching *_lf_0 fluid field is kept too, if present)
mesh_fields = ["integrated_charge", "electrostatic_potential", "E_0"]

# Variance comparisons: [case, estimator], estimator in particle | variance_reduced | low_fidelity.
# Both sides use the same members (same seeds), so the paired-difference bias check is valid.
#[[comparison]]
#name = "vrpic_vs_pic"
#label = "VR-PIC (VR moments, VR E field) vs PIC (standard moments)"
#baseline = ["pic", "particle"]
#test = ["vrpic", "variance_reduced"]

[plot]
steps = [0, 1, 2, 3, 4, 5, 10, 50]
moments = ["number_density", "bulk_velocity_0", "temperature"]
fields = ["integrated_charge", "electrostatic_potential", "E_0"]
    """

    for dt_time_collision_frequency, use_variance_reduction in itertools.product(dt_times_collision_frequencies, [True, False]):
        input_deck_contents = get_input_deck(dt_time_collision_frequency, use_variance_reduction)
        formatted_test_name = format_test_name(dt_time_collision_frequency, use_variance_reduction)
        yaml = f"sod_shock_tube_{formatted_test_name}.yaml"
        with open(yaml, "w") as input_deck:
            input_deck.write(input_deck_contents)

        study_toml_contents += f"""
[cases."{formatted_test_name}"]
deck = "{yaml}"
label = "{formatted_test_name}"
estimator = "{"variance_reduced" if use_variance_reduction else "particle"}"
        """

        # result = subprocess.run([mfpic_executable, "-i", yaml], capture_output=True, text=True)

        # log_filename = f"{output_directory}/execute_{refinement_level:02}.log"
        # with open(log_filename, "w") as execute_log:
        #     execute_log.write(result.stdout)
        #     execute_log.write(result.stderr)

        # verification.check_fluid_energy_positive_and_constant('Total_Fluid_Energy')
        # os.rename('output_lf_0.csv', f"{output_directory}/output_lf_0_{refinement_level:02}.csv")
        # os.rename(yaml, f"{output_directory}/sod_shock_tube_{refinement_level:02}.yaml")

        # result.check_returncode()

    study_toml_filename = "study.toml"
    with open(study_toml_filename, "w") as study_toml:
        study_toml.write(study_toml_contents)

    vrstudy.main(["all", study_toml_filename, "--force"])


def analyze():
    pass


if __name__ == "__main__":
    import sys

    if "run" in sys.argv[1:]:
        mfpic_executable = sys.argv[2]
        run(mfpic_executable)
    else:
        analyze()
