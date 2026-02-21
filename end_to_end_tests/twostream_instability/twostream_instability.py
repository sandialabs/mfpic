# Two cold electron beams streaming in opposite directions interact electrostatically in a periodic, 1D domain.
# Energy transfers from the kinetic energy of the electron beams into electrostatic energy,
# which grows exponentially at an analytically known growth rate.
# The primary mode of the electrostatic waves has wavelength equal to the length of the domain.
# An immobile ion background is added to satisfy the compatibilty condition for the Poisson problem.
# Problem parameters taken from https://doi.org/10.1016/j.cpc.2022.108569

import numpy as np
from scipy.constants import electron_mass, proton_mass, elementary_charge, epsilon_0

number_density = 1.0e16
bulk_speed = 3.2e5
plasma_frequency = np.sqrt(number_density * elementary_charge**2.0 / electron_mass / epsilon_0)
dt = 2.0 * np.pi / plasma_frequency / 120.0
num_time_steps = 450
dx = bulk_speed * dt
num_elements = 400
domain_length = num_elements * dx
num_macroparticles_per_population = 5_000

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
      Number of Macroparticles per Species: {2 * num_macroparticles_per_population}
      Constant:
        Temperature: 0.0
        Number Density: {number_density}
    - Species: [electron]
      Number of Macroparticles per Species: {num_macroparticles_per_population}
      Constant:
        Bulk Velocity: [{bulk_speed}, 0.0, 0.0]
        Temperature: 0.0
        Number Density: {number_density / 2.0}
    - Species: [electron]
      Number of Macroparticles per Species: {num_macroparticles_per_population}
      Constant:
        Bulk Velocity: [-{bulk_speed}, 0.0, 0.0]
        Temperature: 0.0
        Number Density: {number_density / 2.0}

Output:
  Stride: 1

  """
  yaml = "twostream_instability.yaml"
  with open("twostream_instability.yaml", 'w') as input_deck:
    input_deck.write(input_deck_contents)

  result = subprocess.run([mfpic_executable, "-i", yaml])
  result.check_returncode()

def analyze():
  output = np.genfromtxt("output.csv", names=True)

  simulation_times = output["Time"]
  log_energy = np.log10(output["Field_Energy"])

  # Find a stretch of time where the log of the electrostatic energy is "most" linear,
  # judged using the residual for the best fit line.
  # The slope is the growth rate.
  linear_regime_min_timestep_window = 150
  linear_regime_max_timestep_window = 300
  best_linear_regime_fit_residual = 1.0e100
  for timestep_window in range(linear_regime_min_timestep_window, linear_regime_max_timestep_window + 1):
    for timestep_start in range(num_time_steps - timestep_window):
      linear_fit, linear_fit_results = np.polynomial.polynomial.Polynomial.fit(
        simulation_times[timestep_start:timestep_start+timestep_window],
        log_energy[timestep_start:timestep_start+timestep_window],
        1,
        full=True
      )
      residual = linear_fit_results[0][0]
      if residual < best_linear_regime_fit_residual:
        best_linear_regime_timestep_start = timestep_start
        best_linear_regime_timestep_window = timestep_window
        best_linear_regime_fit_residual = residual
        best_linear_regime_fit = linear_fit
  best_fit_growth_rate = best_linear_regime_fit.convert().coef[1]

  wavenumber = 2.0 * np.pi / domain_length
  expected_growth_rate = np.imag(np.sqrt(
    0.5 * (plasma_frequency**2.0 + 2.0 * wavenumber**2.0 * bulk_speed**2.0) -
    0.5 * plasma_frequency * np.sqrt(plasma_frequency**2.0 + 8.0 * wavenumber**2.0 * bulk_speed**2.0) +
    0j
  ))

  assert np.isclose(expected_growth_rate, best_fit_growth_rate, rtol=2.5e-1), f"Expected growth rate of {expected_growth_rate}, but best fit from simulation was {best_fit_growth_rate}"

  return best_linear_regime_timestep_start, best_linear_regime_timestep_window, best_linear_regime_fit, expected_growth_rate

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
