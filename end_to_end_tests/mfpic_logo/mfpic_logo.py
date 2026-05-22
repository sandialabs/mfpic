import numpy as np
from scipy.constants import electron_mass, proton_mass, elementary_charge, epsilon_0

number_density = 1.0e16
plasma_frequency = np.sqrt(number_density * elementary_charge**2.0 / electron_mass / epsilon_0)
dt = 2.0 * np.pi / plasma_frequency / 12000.0
num_time_steps = 450
num_x_cells = 113
num_y_cells = 50
dx = 1.0e-6
num_macroparticles_per_population = 50_000

def run(mfpic_executable):
  import subprocess

  input_deck_contents = f"""
Fields:
  Basis Order: 1

Mesh:
  Type: quad
  Lengths: [{num_x_cells * dx}, {num_y_cells * dx}]
  Number of Elements: [{num_x_cells}, {num_y_cells}]
  Periodic Dimensions: [x, y]

Time Stepping:
  Number of Time Steps: {num_time_steps}
  Time Step Size: {dt}

Species:
  electron:
    Mass: {electron_mass}
    Charge: {-elementary_charge}
  immobile_proton:
    Mass: {proton_mass}
    Charge: {elementary_charge}
    Charge Over Mass: 0.0

Particles:
  Boundary Conditions: []
  Default Boundary Condition: Reflecting
  Initial Conditions:
    - Species: [immobile_proton]
      Number of Macroparticles per Species: {num_macroparticles_per_population}
      Constant:
        Temperature: 0.0
        Number Density: {number_density}
    - Species: [electron]
      Number of Macroparticles per Species: {num_macroparticles_per_population}
      From File:

Output:
  Stride: 1

  """
  yaml = "twostream_instability.yaml"
  with open("twostream_instability.yaml", 'w') as input_deck:
    input_deck.write(input_deck_contents)

  result = subprocess.run([mfpic_executable, "-i", yaml])
  result.check_returncode()

def analyze():
  import matplotlib.pyplot as plt
  from PIL import Image
  import h5py

  x_bins = np.linspace(0.0, dx * num_x_cells, num_x_cells)
  y_bins = np.linspace(0.0, dx * num_y_cells, num_y_cells)

  particle_file = h5py.File("particles.h5part")
  image_file_basename = "electron_number_density_step"
  for timestep in range(num_time_steps + 1):
    particle_timestep_data = particle_file[f"Step#{timestep}"]
    electron_x = particle_timestep_data["x"][num_macroparticles_per_population:]
    electron_y = particle_timestep_data["y"][num_macroparticles_per_population:]
    electron_weights = particle_timestep_data["weight"][num_macroparticles_per_population:]
    fig, ax = plt.subplots()
    ax.hist2d(electron_x, electron_y, bins = [x_bins, y_bins], weights = electron_weights)
    ax.set_axis_off()
    fig.savefig(f"{image_file_basename}_{timestep}.png", dpi=300, bbox_inches="tight")
    plt.close('all')

  reversed_images = [Image.open(f"{image_file_basename}_{frame}.png").convert("RGBA") for frame in reversed(range(num_time_steps + 1))]
  reversed_images = reversed_images[::3]
  first, rest = reversed_images[0], reversed_images[1:]
  duration_ms = 1000.0 / 24.0
  first.save(
    f"{image_file_basename}.gif",
    save_all=True,
    append_images = rest,
    duration = duration_ms,
    loop = 1,
    disposal = 2
  )

if __name__ == "__main__":
  import sys

  if "run" in sys.argv[1:]:
    run(sys.argv[2])
  else:
    analyze()
