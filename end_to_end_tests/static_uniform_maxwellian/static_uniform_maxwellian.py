import sys
sys.path.append("../python")

import numpy as np
import read_mesh_data
from scipy.constants import electron_mass, proton_mass, elementary_charge, epsilon_0, electron_volt, Boltzmann

number_density = 1.0e20
temperature = 10.0 * electron_volt / Boltzmann
debye_length = np.sqrt(epsilon_0 * Boltzmann * temperature / (number_density * elementary_charge**2.0))
dx = debye_length
num_elements = 60
domain_length = num_elements * debye_length
most_probable_speed = np.sqrt(2.0 * Boltzmann * temperature / proton_mass)
acoustic_transit_time = domain_length / most_probable_speed
plasma_frequency = np.sqrt(number_density * elementary_charge**2.0 / proton_mass / epsilon_0)
dt = min(dx / (3.0 * most_probable_speed), 0.1 / plasma_frequency)
num_macroparticles_per_population = 100 * num_elements
num_timesteps = int(np.ceil(acoustic_transit_time / dt))

def run(mfpic_executable):
  import subprocess

  input_deck_contents = f"""
Fields:
  Basis Order: 1
  Boundary Conditions:

Mesh:
  Type: line
  Lengths: [{domain_length}]
  Number of Elements: [{num_elements}]
  Periodic Dimensions: ["x"]

Time Stepping:
  Final Time: {acoustic_transit_time}
  Time Step Size: {dt}

Species:
  proton:
    Mass: {proton_mass}
    Charge: {elementary_charge}

Particles:
  Boundary Conditions:
  Default Boundary Condition: Reflecting
  Initial Conditions:
    - Species: [proton]
      Number of Macroparticles per Species: {num_macroparticles_per_population}
      Constant:
        Temperature: {temperature}
        Number Density: {number_density}

Output:
  Stride: 1

  """

  yaml = "static_uniform_maxwellian.yaml"
  with open("static_uniform_maxwellian.yaml", 'w') as input_deck:
    input_deck.write(input_deck_contents)

  result = subprocess.run([mfpic_executable, "-i", yaml])
  result.check_returncode()

def plot():
  import h5py
  import matplotlib.pyplot as plt

  particle_file = h5py.File("particles.h5part")

  for step in range(num_timesteps):
    particle_data = particle_file[f"/Step#{step}"]

if __name__ == "__main__":
  if "run" in sys.argv[1:]:
    run(sys.argv[2])
  else:
    plot()
