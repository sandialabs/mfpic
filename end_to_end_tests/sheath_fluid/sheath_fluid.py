import sys
sys.path.append("../python")

import numpy as np
import read_mesh_data
import euler, species
from scipy.constants import electron_mass, proton_mass, elementary_charge, epsilon_0, electron_volt, Boltzmann

number_density = 1.0e20
temperature = 10.0 * electron_volt / Boltzmann
debye_length = np.sqrt(epsilon_0 * Boltzmann * temperature / (number_density * elementary_charge**2.0))
dx = debye_length
num_elements = 60
domain_length = num_elements * debye_length
most_probable_ion_speed = np.sqrt(2.0 * Boltzmann * temperature / proton_mass)
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

Euler Fluids:
  Basis Order: 0
  Boundary Conditions:
    - Side: left
      Type: Absorbing
    - Side: right
      Type: Reflecting
  Default Boundary Condition: Reflecting
  Initial Conditions:
    - Species: [proton]
      Constant:
        Temperature: {temperature}
        Number Density: {number_density}
  Sources:
    - Species: [proton]
      Constant:
        Temperature: {temperature}
        Number Density: {source_number_density / dt}

Output:
  Stride: {num_timesteps}

  """

  yaml = "sheath.yaml"
  with open("sheath.yaml", 'w') as input_deck:
    input_deck.write(input_deck_contents)

  result = subprocess.run([mfpic_executable, "-i", yaml])
  result.check_returncode()

def plot():
  import matplotlib.pyplot as plt

  plt.rcParams['text.usetex'] = True

  _, mesh_data = read_mesh_data.read_mesh_data()
  mesh_data_at_last_timestep = mesh_data[-1]

  fluid_data = np.transpose(mesh_data_at_last_timestep["species_0_lf_0"])
  ion_species = species.Species(elementary_charge, proton_mass)

  fig,ax = plt.subplots(1,3,figsize=(12,6),sharex=True)

  temperature = np.zeros_like(fluid_data[0])

  for i in range(0,temperature.shape[0]):
    state = euler.construct_conservative_state(fluid_data[0][i],[fluid_data[1][i],fluid_data[2][i],fluid_data[3][i]],fluid_data[4][i])
    temperature[i] = euler.get_temperature_from_conservative_state(state,ion_species)

  ax[0].plot(mesh_data_at_last_timestep["points"][:,0], fluid_data[0] / ion_species.mass)
  ax[1].plot(mesh_data_at_last_timestep["points"][:,0], fluid_data[1] / fluid_data[0])
  ax[2].plot(mesh_data_at_last_timestep["points"][:,0], temperature)

  ax[0].set_xlabel(r"$x$ $[m]$")
  ax[0].set_ylabel(r"$n$ $[-/m^3]$")
  ax[1].set_ylabel(r"$u_x$ $[m/s]$")
  ax[2].set_ylabel(r"$T$ $[K]$")
  fig.tight_layout()
  fig.savefig("sheath.pdf")

if __name__ == "__main__":
  if "run" in sys.argv[1:]:
    run(sys.argv[2])
  elif "plot" in sys.argv[1:]:
    plot()
