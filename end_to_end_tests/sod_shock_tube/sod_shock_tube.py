from scipy.constants import electron_mass, elementary_charge
import subprocess

domain_length = 1.0
num_elements = 400
dx = domain_length / num_elements
cfl = 0.5
max_wavespeed = 1.
dt = cfl * dx / max_wavespeed
initial_final_time = 0.5 * domain_length / max_wavespeed
# num_time_steps = int(initial_final_time / dt)
# final_time = num_time_steps * dt

num_time_steps = 1
final_time = num_time_steps * dt

input_deck_contents = f"""
Fields:
  Basis Order: 1
  Boundary Conditions: []

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
      Constant:
        Number Density: 1e16
        Temperature: 300

Particles:
  Boundary Conditions: []
  Default Boundary Condition: Reflecting
  Initial Conditions:

Output:
  Stride: 1

"""

def run(mfpic_executable):
    yaml = "sod_shock_tube.yaml"
    with open(yaml, 'w') as input_deck:
        input_deck.write(input_deck_contents)

    result = subprocess.run([mfpic_executable, "-i", yaml])
    result.check_returncode()

def analyze():
    pass

if __name__ == "__main__":
  import sys

  if "run" in sys.argv[1:]:
    run(sys.argv[2])
  else:
    analyze()