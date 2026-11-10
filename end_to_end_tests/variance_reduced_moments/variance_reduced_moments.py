import sys

sys.path.append("../python")

import os
import numpy as np
from scipy.constants import neutron_mass

# A uniform Maxwellian loaded both as particles and as an Euler fluid with the same moments. At the
# initial output f_LF / f_i = 1 for every particle, so the variance reduced number density is the
# control variate integral, i.e. exactly the fluid density, while the standard estimator is noisy.
number_density = 1.0e20
temperature = 300.0
domain_length = 1.0
num_elements = 20
particles_per_element = 100
dt = 1.0e-6
species_name = "neutron"

def run(mfpic_executable):
  import subprocess

  input_deck_contents = f"""
Fields:
  Basis Order: 1

Mesh:
  Type: line
  Lengths: [{domain_length}]
  Number of Elements: [{num_elements}]
  Periodic Dimensions: [x]

Time Stepping:
  Number of Time Steps: 1
  Time Step Size: {dt}

Species:
  {species_name}:
    Mass: {neutron_mass}
    Charge: 0.0

Particles:
  Boundary Conditions: []
  Default Boundary Condition: Reflecting
  Initial Conditions:
    - Species: [{species_name}]
      Number of Macroparticles per Species: {num_elements * particles_per_element}
      Constant:
        Temperature: {temperature}
        Number Density: {number_density}

Euler Fluids:
  Basis Order: 0
  Initial Conditions:
    - Species: [{species_name}]
      Constant:
        Temperature: {temperature}
        Number Density: {number_density}

Variance Reduction:
  Strategy: Euler Fluid

Output:
  Stride: 1
  Particles: False
  Mesh Data: False
  """
  yaml = "input.yaml"
  with open(yaml, 'w') as input_deck:
    input_deck.write(input_deck_contents)

  result = subprocess.run([mfpic_executable, "-i", yaml])
  result.check_returncode()

def read_initial_number_density(prefix):
  filename = f"{prefix}_{species_name}.csv"
  assert os.path.isfile(filename), f"{filename} was not written"
  data = np.genfromtxt(filename, delimiter=",", names=True)
  return data["number_density"][data["step"] == 0]

def analyze():
  standard = read_initial_number_density("particle_moments")
  variance_reduced = read_initial_number_density("variance_reduced_particle_moments")
  fluid = read_initial_number_density("low_fidelity_moments")

  fluid_error = np.max(np.abs(fluid / number_density - 1.0))
  assert fluid_error < 1.0e-12, f"Fluid number density is not the initial condition (relative error {fluid_error})"

  # The standard estimator has ~1/sqrt(particles_per_element) = 10% noise. If the fluid were not used
  # as the control variate, the variance reduced estimate would carry the same noise.
  standard_spread = np.std(standard / number_density)
  assert standard_spread > 1.0e-2, f"Standard estimator unexpectedly noise free (relative std {standard_spread})"

  variance_reduced_error = np.max(np.abs(variance_reduced / fluid - 1.0))
  assert variance_reduced_error < 1.0e-10, f"""
Variance reduced number density does not equal the fluid's for an exact control variate
(max relative difference {variance_reduced_error}; standard estimator relative std {standard_spread}).
  """

if __name__ == "__main__":
  import sys

  if "run" in sys.argv[1:]:
    run(sys.argv[2])
  else:
    analyze()
