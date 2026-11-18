import sys

sys.path.append("../python")

import itertools
import numpy as np
import os
from scipy.constants import electron_mass, proton_mass, elementary_charge, epsilon_0, electron_volt, Boltzmann
import vrstudy

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
particles_per_cell = 100
num_macroparticles_per_population = particles_per_cell * num_elements
bohm_speed = np.sqrt(Boltzmann * temperature / proton_mass)
ion_thermal_speed = bohm_speed
bohm_current_density = number_density * elementary_charge * bohm_speed * np.exp(-0.5)
ion_wall_flux = bohm_current_density / elementary_charge
right_wall = "reservoir"  # "reservoir": absorbing for particles and fluid, with a Maxwellian influx; "reflecting": reflecting particles, influx only
reservoir_outflux_at_initial_density = number_density * ion_thermal_speed / np.sqrt(2.0 * np.pi)
right_wall_influx = ion_wall_flux + (reservoir_outflux_at_initial_density if right_wall == "reservoir" else 0.0)
source_number_density = right_wall_influx / domain_length * dt
source_num_macroparticles = int(source_number_density / number_density * num_macroparticles_per_population)
fluid_volumetric_source_rate = ion_wall_flux / domain_length
fluid_wall_cell_source_rate = ion_wall_flux / dx
fluid_wall_cell_source_bulk_velocity = -ion_thermal_speed * np.sqrt(np.pi / 2.0)
fluid_wall_cell_source_temperature = (2.0 / 3.0) * (2.0 - np.pi / 4.0) * temperature
negligible_source_rate_fraction = 1.0e-12
fluid_source_shape = "boundary"  # "boundary": Maxwellian influx BC matching the particle boundary source; "wall_cell"/"uniform": volumetric sources
particle_source_type = "boundary"  # "volumetric" mirrors the fluid source; "boundary" injects a Maxwellian flux at the right wall
electron_thermal_speed = np.sqrt(8.0 * Boltzmann * temperature / np.pi / electron_mass)
electron_boundary_flux_per_unit_number_density = electron_thermal_speed / 4.0
electron_source_rate = ion_wall_flux
electron_reference_number_density = electron_source_rate / electron_boundary_flux_per_unit_number_density
num_timesteps = int(np.ceil(ion_acoustic_transit_time / dt))
output_stride = 5

dt_times_collision_frequencies = [0.1]
pdf_update_stride = None  # set to an int to enable periodic f_HF re-estimation


def get_volumetric_source_state(number_density_scale):
    if fluid_source_shape == "uniform":
        return f"""      Constant:
        Temperature: {temperature}
        Number Density: {number_density_scale * fluid_volumetric_source_rate}"""
    if fluid_source_shape == "wall_cell":
        return f"""      Sod:
        Discontinuity Location: {domain_length - dx}
        Left State:
          Temperature: {temperature}
          Number Density: {number_density_scale * negligible_source_rate_fraction * fluid_wall_cell_source_rate}
          Bulk Velocity: [0.0, 0.0, 0.0]
        Right State:
          Temperature: {fluid_wall_cell_source_temperature}
          Number Density: {number_density_scale * fluid_wall_cell_source_rate}
          Bulk Velocity: [{fluid_wall_cell_source_bulk_velocity}, 0.0, 0.0]"""
    raise ValueError(f"Unknown fluid source shape: {fluid_source_shape}")


def get_fluid_source_block():
    if fluid_source_shape == "boundary":
        return ""
    return f"""  Sources:
    - Species: [proton]
{get_volumetric_source_state(1.0)}"""


def get_particle_right_boundary_block():
    if right_wall == "reservoir":
        return """    - Side: right
      Type: Absorbing
"""
    return ""


def get_fluid_right_boundary_block():
    if fluid_source_shape == "boundary":
        return f"""    - Side: right
      Type: Maxwellian Influx
      Number Flux Density: {right_wall_influx}
      Temperature: {temperature}"""
    return """    - Side: right
      Type: Reflecting"""


def get_particle_source_block():
    if particle_source_type == "volumetric":
        return f"""  Sources:
    - Species: [proton]
      Number of Macroparticles per Species: {source_num_macroparticles}
{get_volumetric_source_state(dt)}"""
    if particle_source_type == "boundary":
        return f"""  Boundary Sources:
    - Species: [proton]
      Side: right
      Number of Macroparticles per Species per Timestep: {source_num_macroparticles}
      Maxwellian:
        Temperature: {temperature}
        Flux Density: {right_wall_influx}"""
    raise ValueError(f"Unknown particle source type: {particle_source_type}")


def format_test_name(dt_times_collision_frequency: float, use_variance_reduction: bool):
    return f'nu_x_dt_{dt_times_collision_frequency}_vr_{use_variance_reduction}'


def get_input_deck(dt_times_collision_frequency: float, use_variance_reduction: bool):
    collision_frequency = dt_times_collision_frequency / dt

    no_variance_reduction = """
Variance Reduction:
  Strategy: none
    """

    with_variance_reduction = """
Variance Reduction:
  Strategy: Euler Fluid
  Cell Limiting: Every Particle
    """

    pdf_update_block = "" if pdf_update_stride is None else f"""
PDF Update:
  Stride: {pdf_update_stride}
"""

    collision_block = f"""
Collisions:
    - BGK Relaxation:
        Species: proton
        Collision Frequency: {collision_frequency}
    """

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
  Number of Time Steps: {num_timesteps}
  Time Step Size: {dt}

Species:
  proton:
    Mass: {proton_mass}
    Charge: {elementary_charge}

Particles:
  Boundary Conditions:
    - Side: left
      Type: Absorbing
{get_particle_right_boundary_block()}  Default Boundary Condition: Reflecting
  Initial Conditions:
    - Species: [proton]
      Number of Macroparticles per Species: {num_macroparticles_per_population}
      Constant:
        Temperature: {temperature}
        Number Density: {number_density}
{get_particle_source_block()}

Euler Fluids:
  Basis Order: 0
  Initial Conditions:
    - Species: [proton]
      Constant:
        Temperature: {temperature}
        Number Density: {number_density}
{get_fluid_source_block()}
  Boundary Conditions:
    - Side: left
      Type: Absorbing
{get_fluid_right_boundary_block()}

{collision_block if collision_frequency > 0 else ""}

{pdf_update_block}

Output:
  Stride: {output_stride}
  Mesh Output Folder: MeshOutput_{format_test_name(dt_times_collision_frequency, use_variance_reduction)}

{with_variance_reduction if use_variance_reduction else no_variance_reduction}

"""
    return input_deck_contents


def run(mfpic_executable):
    output_directory = "./"
    os.makedirs(output_directory, exist_ok=True)

    plot_steps = sorted({0, output_stride, 2 * output_stride, 10 * output_stride, 20 * output_stride,
                         40 * output_stride, (num_timesteps // output_stride) * output_stride})

    study_toml_contents = f"""
[study]
name = "proton"
mfpic = "{mfpic_executable}"              # executable name on PATH, or an absolute path
members = 112                # ensemble size per case
base_seed = 1000             # member i uses Particles: Seed = base_seed + i in every case
jobs = 112                   # concurrent mfpic runs (each pinned to one core, OMP_NUM_THREADS=1)
pin_cores = true
species = ["proton"]
keep_particle_dumps = [1]    # members that keep particles.h5part (all others turn it off)
keep_mesh_output = true

[collect]
# MeshOutput point fields to keep (the matching *_lf_0 fluid field is kept too, if present)
mesh_fields = ["integrated_charge", "electrostatic_potential", "E_0"]

[plot]
steps = {plot_steps}
moments = ["number_density", "bulk_velocity_0", "temperature"]
fields = ["integrated_charge", "electrostatic_potential", "E_0"]
    """

    for dt_times_collision_frequency, use_variance_reduction in itertools.product(dt_times_collision_frequencies, [True, False]):
        input_deck_contents = get_input_deck(dt_times_collision_frequency, use_variance_reduction)
        formatted_test_name = format_test_name(dt_times_collision_frequency, use_variance_reduction)
        yaml = f"sheath_{formatted_test_name}.yaml"
        with open(yaml, "w") as input_deck:
            input_deck.write(input_deck_contents)

        study_toml_contents += f"""
[cases."{formatted_test_name}"]
deck = "{yaml}"
label = "{formatted_test_name}"
estimator = "{"variance_reduced" if use_variance_reduction else "particle"}"
        """

    for dt_times_collision_frequency in dt_times_collision_frequencies:
        vr_case = format_test_name(dt_times_collision_frequency, True)
        pic_case = format_test_name(dt_times_collision_frequency, False)
        study_toml_contents += f"""
[[comparison]]
name = "vrpic_vs_pic_{dt_times_collision_frequency}"
label = "VR-PIC (VR moments) vs PIC (standard moments), nu dt = {dt_times_collision_frequency}"
baseline = ["{pic_case}", "particle"]
test = ["{vr_case}", "variance_reduced"]

[[comparison]]
name = "vr_estimator_{dt_times_collision_frequency}"
label = "VR moments vs standard moments of the same VR-PIC run, nu dt = {dt_times_collision_frequency}"
baseline = ["{vr_case}", "particle"]
test = ["{vr_case}", "variance_reduced"]
        """

    study_toml_filename = "study.toml"
    with open(study_toml_filename, "w") as study_toml:
        study_toml.write(study_toml_contents)

    vrstudy.main(["all", study_toml_filename, "--force"])


def analyze():
    pass


if __name__ == "__main__":
    if "run" in sys.argv[1:]:
        mfpic_executable = sys.argv[2]
        run(mfpic_executable)
    else:
        analyze()
