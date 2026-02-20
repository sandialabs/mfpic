import euler
import matplotlib.pyplot as plt
import numpy as np
import os
import read_mesh_data
from scipy.constants import electron_mass
import species
import subprocess

electron_species = species.Species(mass = electron_mass)

mass_density_left = 1.0
pressure_left = 1.0
bulk_velocity_left = np.array([0., 0., 0.])

mass_density_right = 0.125
pressure_right = 0.1
bulk_velocity_right = np.array([0., 0., 0.])

number_density_left = mass_density_left / electron_species.mass
temperature_left = euler.temperature(number_density_left, pressure_left)

number_density_right = mass_density_right / electron_species.mass
temperature_right = euler.temperature(number_density_right, pressure_right)

sound_speed_left = euler.speed_of_sound(electron_species, mass_density_left, pressure_left)
max_wavespeed_left = np.max(np.abs(bulk_velocity_left)) + sound_speed_left

sound_speed_right = euler.speed_of_sound(electron_species, mass_density_right, pressure_right)
max_wavespeed_right = np.max(np.abs(bulk_velocity_right)) + sound_speed_right

max_wavespeed = max(max_wavespeed_left, max_wavespeed_right)

domain_length = 1.0
num_elements = 400
dx = domain_length / num_elements
cfl = 0.5
dt = cfl * dx / max_wavespeed
initial_final_time = 0.5 * domain_length / max_wavespeed
num_time_steps = int(initial_final_time / dt)
final_time = num_time_steps * dt

input_deck_contents = f"""
Fields:
  Basis Order: 1
  Boundary Conditions: []

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
        Discontinuity Location: {0.5 * domain_length}
        Left State:
          Number Density: {number_density_left}
          Temperature: {temperature_left}
          Bulk Velocity: [{bulk_velocity_left[0]}, {bulk_velocity_left[1]}, {bulk_velocity_left[2]}]
        Right State:
          Number Density: {number_density_right}
          Temperature: {temperature_right}
          Bulk Velocity: [{bulk_velocity_right[0]}, {bulk_velocity_right[1]}, {bulk_velocity_right[2]}]

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
  timesteps, mesh_data = read_mesh_data.read_mesh_data()

  points = mesh_data[0]['points']
  x_points = points[:, 0]

  os.makedirs("Figures", exist_ok=True)
  for i in range(len(timesteps)):
    mass_density = mesh_data[i]['species_0'][:, 0]
    time = timesteps[i]

    fig, axes = plt.subplots()
    axes.plot(x_points, mass_density)
    axes.set_title(f"Mass Density At Time = {time}")
    axes.set_xlabel('x')
    axes.set_ylabel('rho')
    fig.savefig(f"Figures/MassDensity{i:03}.png")
    plt.close(fig)

if __name__ == "__main__":
  import sys

  if "run" in sys.argv[1:]:
    run(sys.argv[2])
  else:
    analyze()