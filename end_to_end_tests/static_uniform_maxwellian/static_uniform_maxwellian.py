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
final_time = 2*acoustic_transit_time
plasma_frequency = np.sqrt(number_density * elementary_charge**2.0 / proton_mass / epsilon_0)
dt = min(dx / (3.0 * most_probable_speed), 0.1 / plasma_frequency)
num_macroparticles_per_population = 100 * num_elements
num_timesteps = int(np.ceil(final_time / dt))

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
  Final Time: {final_time}
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
  import matplotlib.animation as animation
  import matplotlib.pyplot as plt

  particle_file = h5py.File("particles.h5part")

  particle_data = particle_file[f"/Step#0"]
  particle_skip_interval = 20
  velocities = particle_data["vx"][::particle_skip_interval]
  pdfs = particle_data["particle_distribution_function_value"][::particle_skip_interval]

  fig, ax = plt.subplots()
  xlim = (-2.5, 2.5)
  ax.set_xlim(xlim)
  ax.set_ylabel(r"$f$")
  ax.set_xlabel(r"$v_x / v_{th}$")
  phase_space_scatter = ax.scatter(velocities / most_probable_speed, pdfs / number_density, color="C1")
  plot_points = np.linspace(xlim[0], xlim[1], 100)
  ax.plot(plot_points, np.sqrt(proton_mass / (2.0 * np.pi * Boltzmann * temperature)) * np.exp(- np.power(plot_points, 2.0)), color="C2")

  def update_phase_space(frame):
    particle_data = particle_file[f"/Step#{frame}"]
    velocities = particle_data["vx"][::particle_skip_interval]
    pdfs = particle_data["particle_distribution_function_value"][::particle_skip_interval]
    updated_phase_space_data = np.stack([velocities / most_probable_speed, pdfs / number_density]).T
    phase_space_scatter.set_offsets(updated_phase_space_data)
    return phase_space_scatter

  phase_space_animation = animation.FuncAnimation(
    fig=fig,
    func=update_phase_space,
    frames=num_timesteps+1,
    interval=30
  )
  phase_space_animation.save("phase_space.mp4")

if __name__ == "__main__":
  if "run" in sys.argv[1:]:
    run(sys.argv[2])
  else:
    plot()
