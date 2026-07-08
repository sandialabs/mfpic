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
most_probable_ion_speed = np.sqrt(2.0 * Boltzmann * temperature / proton_mass)
ion_thermal_speed = np.sqrt(Boltzmann * temperature / proton_mass)
ion_acoustic_transit_time = domain_length / most_probable_ion_speed
ion_plasma_frequency = np.sqrt(number_density * elementary_charge**2.0 / proton_mass / epsilon_0)
dt = min(dx / (3.0 * most_probable_ion_speed), 0.1 / ion_plasma_frequency)
num_macroparticles_per_population = 100 * num_elements
bohm_speed = np.sqrt(Boltzmann * temperature / proton_mass)
bohm_current_density = number_density * elementary_charge * bohm_speed * np.exp(-0.5)
ion_wall_flux = bohm_current_density / elementary_charge
source_number_density = ion_wall_flux / domain_length * dt
source_num_macroparticles = int(source_number_density / number_density * num_macroparticles_per_population)
electron_thermal_speed = np.sqrt(8.0 * Boltzmann * temperature / np.pi / electron_mass)
electron_boundary_flux_per_unit_number_density = electron_thermal_speed / 4.0
electron_source_rate = ion_wall_flux
electron_reference_number_density = electron_source_rate / electron_boundary_flux_per_unit_number_density
num_timesteps = int(np.ceil(ion_acoustic_transit_time / dt))

def run(mfpic_executable):
  import subprocess

  input_deck_contents = f"""
Fields:
  Basis Order: 1
  Boltzmann Electrons:
    Reference Number Density: {electron_reference_number_density}
    Temperature: {temperature}
  Boundary Conditions:
    - Side: left
      Value: 0.0

Mesh:
  Type: line
  Lengths: [{domain_length}]
  Number of Elements: [{num_elements}]

Time Stepping:
  Final Time: {ion_acoustic_transit_time}
  Time Step Size: {dt}

Species:
  proton:
    Mass: {proton_mass}
    Charge: {elementary_charge}

Particles:
  Boundary Conditions:
    - Side: left
      Type: Absorbing
  Default Boundary Condition: Reflecting
  Initial Conditions:
    - Species: [proton]
      Number of Macroparticles per Species: {num_macroparticles_per_population}
      Constant:
        Temperature: {temperature}
        Number Density: {number_density}
  Sources:
    - Species: [proton]
      Number of Macroparticles per Species: {source_num_macroparticles}
      Constant:
        Temperature: {temperature}
        Number Density: {source_number_density}

Output:
  Stride: 1

  """

  yaml = "sheath.yaml"
  with open("sheath.yaml", 'w') as input_deck:
    input_deck.write(input_deck_contents)

  result = subprocess.run([mfpic_executable, "-i", yaml])
  result.check_returncode()

def driftVelocityInElement(particle_data, element):
  weights = particle_data["weight"][:]
  vx = particle_data["vx"][:]
  elements = particle_data["element"][:]

  num_physical_particles = np.sum(weights, where = elements == element)
  weighted_velocity = np.sum(weights * vx, where = elements == element)

  return weighted_velocity / num_physical_particles

def findElementWhereIonDriftSpeedSatisfiesBohmCriterion(particle_data):
  for element in range(num_elements):
    if driftVelocityInElement(particle_data, element) > -bohm_speed:
      return element

def analyze():
  import h5py

  particle_file = h5py.File("particles.h5part")
  particle_data_at_last_timestep = particle_file[f"Step#{num_timesteps}"]
  sheath_entrance_element = findElementWhereIonDriftSpeedSatisfiesBohmCriterion(particle_data_at_last_timestep)
  sheath_entrance_node = sheath_entrance_element - 1

  _, mesh_data = read_mesh_data.read_mesh_data()
  mesh_data_at_last_timestep = mesh_data[-1]
  simulated_sheath_entrance_potential = mesh_data_at_last_timestep["electrostatic_potential"][2*sheath_entrance_node]

  expected_sheath_entrance_potential = Boltzmann * temperature / (2.0 * elementary_charge) * np.log(
    proton_mass / electron_mass / 4.0 / np.pi
  )

  assert np.isclose(expected_sheath_entrance_potential, simulated_sheath_entrance_potential, rtol=2.5e-1), f"Expected potential at sheath entrance to be {expected_sheath_entrance_potential} V, but computed {simulated_sheath_entrance_potential} V"

def plot():
  import h5py
  import matplotlib.animation as animation
  import matplotlib.pyplot as plt
  import matplotlib as mpl

  particle_moments = np.genfromtxt("particle_moments_proton.csv", delimiter=",", names=True)
  cell_centers = np.compress(particle_moments["step"] == 0, particle_moments["x"])
  number_densities = np.compress(particle_moments["step"] == 0, particle_moments["number_density"])
  number_density_fig, number_density_ax = plt.subplots()
  number_density_plot, = number_density_ax.plot(
    cell_centers / debye_length,
    number_densities
  )
  number_density_ax.set_xlim((cell_centers[0] / debye_length, cell_centers[-1] / debye_length))
  number_density_ax.set_xlabel(r"$x / \lambda_D$")
  number_density_ax.set_ylim((0, number_density * 1.4))
  number_density_ax.set_ylabel(r"Number density ($m^{-3}$)")

  def update_number_density(frame):
    number_densities = np.compress(particle_moments["step"] == frame, particle_moments["number_density"])
    number_density_plot.set_ydata(number_densities)
    return number_density_plot

  number_density_animation = animation.FuncAnimation(
    fig=number_density_fig,
    func=update_number_density,
    frames=num_timesteps,
    interval=30,
  )
  number_density_animation.save("number_density.mp4")

  particle_file = h5py.File("particles.h5part")

  _, mesh_data = read_mesh_data.read_mesh_data()
  mesh_data_at_last_timestep = mesh_data[-1]
  mesh_points = mesh_data_at_last_timestep["points"]

  velocity_colormap = mpl.colormaps["plasma"]
  velocity_min = -4
  velocity_max =  4
  velocity_colormap_normalizer = mpl.colors.Normalize(vmin=velocity_min, vmax=velocity_max)

  phase_space_fig, phase_space_ax = plt.subplots()
  particle_data = particle_file["/Step#0"]
  phase_space_scatter = phase_space_ax.scatter(
    particle_data["x"] / debye_length,
    particle_data["vx"] / ion_thermal_speed,
    marker='.',
    c=particle_data["vx"] / ion_thermal_speed,
    cmap=velocity_colormap,
    norm=velocity_colormap_normalizer
  )
  phase_space_ax.set_xlim((mesh_points[0][0] / debye_length, mesh_points[-1][0] / debye_length))
  phase_space_ax.set_xlabel(r"$x / \lambda_D$")
  phase_space_ax.set_ylim((velocity_min, velocity_max))
  phase_space_ax.set_ylabel(r"$v_x / v_{th}$")

  def update_phase_space(frame):
    particle_data = particle_file[f"/Step#{frame}"]
    updated_phase_space_data = np.stack([particle_data["x"] / debye_length, particle_data["vx"] / ion_thermal_speed]).T
    phase_space_scatter.set_offsets(updated_phase_space_data)
    phase_space_scatter.set_array(particle_data["vx"] / ion_thermal_speed)
    return phase_space_scatter

  phase_space_animation = animation.FuncAnimation(
    fig=phase_space_fig,
    func=update_phase_space,
    frames=num_timesteps,
    interval=30
  )
  phase_space_animation.save("phase_space.mp4")

if __name__ == "__main__":
  if "run" in sys.argv[1:]:
    run(sys.argv[2])
  elif "plot" in sys.argv[1:]:
    plot()
  else:
    analyze()
