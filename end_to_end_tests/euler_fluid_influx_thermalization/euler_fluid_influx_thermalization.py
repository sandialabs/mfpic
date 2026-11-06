import sys

sys.path.append("../python")

import euler
import numpy as np
import read_mesh_data
from scipy.constants import neutron_mass, Boltzmann
import species
import verification

reservoir_number_density = 1.0e26
reservoir_temperature = 300.0
reservoir_mean_speed = np.sqrt(8.0 * Boltzmann * reservoir_temperature / (np.pi * neutron_mass))
influx_density = reservoir_number_density * reservoir_mean_speed / 4.0
initial_number_density = 0.5 * reservoir_number_density
initial_temperature = 5.0 * reservoir_temperature
domain_length = 1.0
num_elements = 50
dx = domain_length / num_elements
dt = 0.5 * dx / reservoir_mean_speed
transit_time = domain_length / reservoir_mean_speed
simulation_time = 50.0 * transit_time
num_time_steps = int(np.ceil(simulation_time / dt))
dt = simulation_time / num_time_steps
mesh_folder_name = "Output"

def run(mfpic_executable):
  import subprocess

  input_deck_contents = f"""
Fields:
  Basis Order: 1

Mesh:
  Type: line
  Lengths: [{domain_length}]
  Number of Elements: [{num_elements}]

Time Stepping:
  Number of Time Steps: {num_time_steps}
  Time Step Size: {dt}

Species:
  neutron:
    Mass: {neutron_mass}
    Charge: 0.0

Euler Fluids:
  Basis Order: 0
  Initial Conditions:
    - Species: [neutron]
      Constant:
        Temperature: {initial_temperature}
        Number Density: {initial_number_density}
  Boundary Conditions:
    - Side: left
      Type: Maxwellian Influx
      Number Flux Density: {influx_density}
      Temperature: {reservoir_temperature}
    - Side: right
      Type: Maxwellian Influx
      Number Flux Density: {influx_density}
      Temperature: {reservoir_temperature}

Output:
  Stride: {num_time_steps}
  Mesh Output Folder: {mesh_folder_name}

  """
  yaml = "input.yaml"
  with open("input.yaml", 'w') as input_deck:
    input_deck.write(input_deck_contents)

  result = subprocess.run([mfpic_executable, "-i", yaml])
  result.check_returncode()

def analyze():
  neutron_species = species.Species(charge=0.0, mass=neutron_mass)

  _, mesh_data = read_mesh_data.read_mesh_data(mesh_folder_name)
  points = mesh_data[-1]["points"]
  x_points = points[:, 0]

  fluid_data = np.transpose(mesh_data[-1]["species_0_lf_0"])

  relative_tolerance = 0.05

  number_density = euler.get_number_density_from_conservative_state(fluid_data, neutron_species)
  number_density_error = verification.compute_L1_relative_error_1D(
    verification.create_1D_interpolater(number_density, x_points),
    lambda x: reservoir_number_density,
    x_points
  )
  assert(
    number_density_error <= relative_tolerance
  ), f"""
Relative error norm of number density ({number_density_error}) exceeeded tolerance ({relative_tolerance}).
Exact: {reservoir_number_density}
Computed: {number_density}
  """

  temperature = euler.get_temperature_from_conservative_state(fluid_data, neutron_species)
  temperature_error = verification.compute_L1_relative_error_1D(
    verification.create_1D_interpolater(temperature, x_points),
    lambda x: reservoir_temperature,
    x_points
  )
  assert(
    temperature_error <= relative_tolerance
  ), f"""
Relative error norm of temperature ({temperature_error}) exceeeded tolerance ({relative_tolerance}).
Exact: {reservoir_temperature}
Computed: {temperature}
  """

if __name__ == "__main__":
  import sys

  if "run" in sys.argv[1:]:
    run(sys.argv[2])
  else:
    analyze()
