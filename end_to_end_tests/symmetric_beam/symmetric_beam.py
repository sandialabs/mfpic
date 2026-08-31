import sys

sys.path.append("../python")

import species
import utils

import numpy as np
from scipy.constants import (
    Boltzmann,
    electron_mass,
    proton_mass,
    elementary_charge,
    epsilon_0,
)
import subprocess

# ============================================================
# Domain
# ============================================================

domain_length = 0.01
num_elements = 300
dx = domain_length / num_elements

# ============================================================
# Species
# ============================================================

electron_species = species.Species(
    charge=-elementary_charge,
    mass=electron_mass,
)

ion_species = species.Species(
    charge=elementary_charge,
    mass=proton_mass,
)

# ============================================================
# Plasma parameters
# ============================================================

background_number_density = 9.9e15

beam_number_density = 1.0e15

background_temperature = 1000.0

beam_temperature = 500.0

vth = np.sqrt(
    Boltzmann * background_temperature / electron_mass
)

beam_velocity = 4.0 * vth

print(f"vth                    = {vth:.3e} m/s")
print(f"beam velocity          = {beam_velocity:.3e} m/s")
print(
    f"total beam fraction    = "
    f"{2.0 * beam_number_density / background_number_density:.3f}"
)

# ============================================================
# Plasma frequency
# ============================================================

total_electron_density = (
    background_number_density
    + 2.0 * beam_number_density
)

omega_pe = np.sqrt(
    total_electron_density
    * elementary_charge**2
    / (epsilon_0 * electron_mass)
)

print(f"omega_pe               = {omega_pe:.3e} rad/s")

# ============================================================
# Time stepping
# ============================================================

final_time = 100.0 / omega_pe

max_cfl = 0.5

max_wavespeed = (
    beam_velocity
    + 5.0 * vth
)

dt, num_time_steps = (
    utils.compute_timestepping_that_satisfies_cfl(
        max_cfl,
        dx,
        max_wavespeed,
        final_time,
    )
)

print(f"dt                     = {dt:.3e} s")
print(f"num_time_steps         = {num_time_steps}")

# ============================================================
# Macroparticles
# ============================================================

basis_order = 0

num_background_macroparticles = 38000
num_beam_macroparticles = 1000
num_ion_macroparticles = 40000

# ============================================================
# Input deck
# ============================================================

input_deck_contents = f"""
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
    Mass: {electron_species.mass}
    Charge: {electron_species.charge}

Particles:
  Basis Order: {basis_order}

  Initial Conditions:

    # ------------------------------------------------
    # Maxwellian electron background
    # ------------------------------------------------

    - Species: [electron]
      Number of Macroparticles per Species: {num_background_macroparticles}
      Constant:
        Number Density: {background_number_density}
        Temperature: {background_temperature}
        Bulk Velocity: [0.0, 0.0, 0.0]

    # ------------------------------------------------
    # Positive beam
    # ------------------------------------------------

    - Species: [electron]
      Number of Macroparticles per Species: {num_beam_macroparticles}
      Constant:
        Number Density: {beam_number_density}
        Temperature: {beam_temperature}
        Bulk Velocity: [{beam_velocity}, 0.0, 0.0]

    # ------------------------------------------------
    # Negative beam
    # ------------------------------------------------

    - Species: [electron]
      Number of Macroparticles per Species: {num_beam_macroparticles}
      Constant:
        Number Density: {beam_number_density}
        Temperature: {beam_temperature}
        Bulk Velocity: [{-beam_velocity}, 0.0, 0.0]

Output:
  Stride: 50
"""

# ============================================================
# Run
# ============================================================

def run(mfpic_executable):
    yaml = "symmetric_beam.yaml"

    with open(yaml, "w") as input_deck:
        input_deck.write(input_deck_contents)

    result = subprocess.run(
        [mfpic_executable, "-i", yaml]
    )

    result.check_returncode()


def analyze():
    pass


if __name__ == "__main__":
    if "run" in sys.argv[1:]:
        run(sys.argv[2])
    else:
        analyze()
