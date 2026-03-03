import sys

sys.path.append("../python")

import euler
import species
import utils

from scipy.constants import electron_mass
import subprocess

domain_length = 2.0
base_num_elements = 500

gaussian_center = 0.25 * domain_length
gaussian_standard_deviation = 0.01

number_density_offset = 1e19
perturbation = 0.1
number_density_height = perturbation * number_density_offset

temperature = 300
velocity = 3.0
max_cfl = 0.9

final_time = 0.5 * domain_length / velocity

neutral_species = species.Species(mass = electron_mass)
max_mass_density = (number_density_height + number_density_offset) * species.mass
pressure = euler.pressure_from_number_density(number)
sound_speed = euler.speed_of_sound()


def format_mesh_folder_name(refinement_level):
    return f"MeshOutput{refinement_level:02}"

def get_input_deck(refinement_level):
    num_elements = base_num_elements * refinement_level
    dx = domain_length / num_elements
    dt, num_time_steps = utils.compute_timestepping_that_satisfies_cfl(max_cfl, dx, velocity, final_time)

    num_time_steps = 1

    mesh_folder_name = format_mesh_folder_name(refinement_level)

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
  neutral_electron:
    Mass: {electron_mass}
    Charge: 0

Euler Fluids:
  Basis Order: 0
  Initial Conditions:
    - Species: [neutral_electron]
      Gaussian:
        Center: [{gaussian_center}]
        Standard Deviation: {gaussian_standard_deviation}
        Offsets:
          Number Density: {number_density_offset}
          Temperature: {temperature}
          Bulk Velocity: [{velocity}, 0., 0.]
        Heights:
          Number Density: {number_density_height}
          Temperature: 0.

Output:
  Stride: 1
  Mesh Output Folder: {mesh_folder_name}
"""
    return input_deck_contents


def run(mfpic_executable):
    input_deck_contents = get_input_deck(1)
    yaml = "advecting_gaussian.yaml"
    with open(yaml, "w") as input_deck:
        input_deck.write(input_deck_contents)

    print(mfpic_executable, "-i", yaml)
    result = subprocess.run([mfpic_executable, "-i", yaml])
    result.check_returncode()


def analyze():
    pass


def plot():
    pass


if __name__ == "__main__":
    import sys

    if "run" in sys.argv[1:]:
        run(sys.argv[2])
    elif "plot" in sys.argv[1:]:
        plot()
    else:
        analyze()
