import numpy as np
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
num_macroparticles_per_population = 1000 * num_elements
bohm_speed = np.sqrt(Boltzmann * temperature / proton_mass)
bohm_current_density = number_density * elementary_charge * bohm_speed * np.exp(-0.5)
ion_wall_flux = bohm_current_density / elementary_charge
source_number_density = ion_wall_flux / domain_length * dt
source_num_macroparticles = int(source_number_density / number_density * num_macroparticles_per_population)
electron_thermal_speed = np.sqrt(8.0 * Boltzmann * temperature / np.pi / electron_mass)
electron_boundary_flux_per_unit_number_density = electron_thermal_speed / 4.0
electron_source_rate = ion_wall_flux
electron_reference_number_density = electron_source_rate / electron_boundary_flux_per_unit_number_density

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

  """

  yaml = "sheath.yaml"
  with open("sheath.yaml", 'w') as input_deck:
    input_deck.write(input_deck_contents)

  result = subprocess.run([mfpic_executable, "-i", yaml])
  result.check_returncode()

def analyze():
  return

def plot():
  import matplotlib.pyplot as plt

  best_linear_regime_timestep_start, best_linear_regime_timestep_window, best_linear_regime_fit, expected_growth_rate = analyze()

  output = np.genfromtxt("output.csv", names=True)

  simulation_times = output["Time"]
  energy = output["Field_Energy"]

  windowed_times = simulation_times[best_linear_regime_timestep_start:best_linear_regime_timestep_start+best_linear_regime_timestep_window]

  plt.semilogy(simulation_times *plasma_frequency / (2.0 * np.pi), energy, label="simulation")
  plt.semilogy(
    windowed_times *plasma_frequency / (2.0 * np.pi),
    np.pow(10, best_linear_regime_fit(windowed_times)),
    label="growth fit"
  )
  plt.semilogy(
    windowed_times *plasma_frequency / (2.0 * np.pi),
    energy[best_linear_regime_timestep_start] * np.pow(10, expected_growth_rate * (windowed_times - windowed_times[0])),
    label="expected growth"
  )
  plt.axvline(simulation_times[best_linear_regime_timestep_start] *plasma_frequency / (2.0 * np.pi))
  plt.axvline(simulation_times[best_linear_regime_timestep_start+best_linear_regime_timestep_window] *plasma_frequency / (2.0 * np.pi))
  plt.xlabel("Simulation time (plasma frequencies)")
  plt.ylabel("Electrostatic energy (J)")
  plt.legend()
  plt.savefig("twostream.png")

if __name__ == "__main__":
  import sys

  if "run" in sys.argv[1:]:
    run(sys.argv[2])
  elif "plot" in sys.argv[1:]:
    plot()
  else:
    analyze()
