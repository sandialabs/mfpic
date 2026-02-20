import numpy as np
from scipy.constants import Boltzmann as boltzmann_constant

def kinetic_energy_density(mass_density, momentum_density):
    return 0.5 / mass_density * np.dot(momentum_density, momentum_density)

def pressure(number_density, temperature):
    return number_density * temperature * boltzmann_constant

def pressure(species, internal_energy_density):
    return (species.specific_heat_ratio - 1.) * internal_energy_density

def internal_energy_per_unit_mass(species, mass_density, pressure):
    return pressure / ((species.specific_heat_ratio - 1.) * mass_density)

def internal_energy_density(internal_energy_density_per_unit_mass, mass_density):
    return internal_energy_density_per_unit_mass * mass_density

def temperature(number_density, pressure):
    return pressure / (number_density * boltzmann_constant)

def speed_of_sound(species, mass_density, pressure):
    return np.sqrt(species.specific_heat_ratio * pressure / mass_density)

def construct_conservative_state(mass_density, momentum_density, total_energy_density):
    return np.array([mass_density, momentum_density[0], momentum_density[1], momentum_density[2], total_energy_density])

def construct_primitive_state(number_density, bulk_velocity, temperature):
    return np.array([number_density, bulk_velocity[0], bulk_velocity[1], bulk_velocity[2], temperature])

def get_number_density_from_conservative_state(conservative_state, species):
    return conservative_state[0] / species.mass

def get_momentum_density_from_conservative_state(conservative_state):
    return np.array(conservative_state[1:4])

def get_bulk_velocity_from_conservative_state(conservative_state):
    momentum_density = get_momentum_density_from_conservative_state(conservative_state)
    mass_density = conservative_state[0]
    return momentum_density / mass_density

def get_kinetic_energy_density_from_conservative_state(conservative_state):
    mass_density = conservative_state[0]
    momentum_density = get_momentum_density_from_conservative_state(conservative_state)
    return kinetic_energy_density(mass_density, momentum_density)

def get_pressure_from_conservative_state(conservative_state, species):
    total_energy_density = conservative_state[4]
    kinetic_energy_density = get_kinetic_energy_density_from_conservative_state(conservative_state)
    internal_energy_density = total_energy_density - kinetic_energy_density
    return pressure(species, internal_energy_density)

def get_temperature_from_conservative_state(conservative_state, species):
    number_density = get_number_density_from_conservative_state(conservative_state, species)
    pressure = get_pressure_from_conservative_state(conservative_state, species)
    return temperature(number_density, pressure)

def get_mass_density_from_primitive_state(primitive_state, species):
    return primitive_state[0] * species.mass

def get_bulk_velocity_from_primitive_state(primitive_state):
    return np.array(primitive_state[1:4])

def get_pressure_from_primitive_state(primitive_state):
    number_density = primitive_state[0]
    temperature = primitive_state[4]
    return pressure(number_density, temperature)

def get_internal_energy_density_from_primitive_state(primitive_state, species):
    mass_density = get_mass_density_from_primitive_state(primitive_state, species)
    pressure = get_pressure_from_primitive_state(primitive_state)
    internal_energy_per_unit_mass_val = internal_energy_per_unit_mass(species, mass_density, pressure)
    return internal_energy_density(internal_energy_per_unit_mass_val, mass_density)

def convert_from_conservative_to_primitive(conservative_state, species):
    number_density = get_number_density_from_conservative_state(conservative_state, species)
    bulk_velocity = get_bulk_velocity_from_conservative_state(conservative_state)
    temperature = get_temperature_from_conservative_state(conservative_state, species)
    return construct_primitive_state(number_density, bulk_velocity, temperature)

def convert_from_primitive_to_conservative(primitive_state, species):
    mass_density = get_mass_density_from_primitive_state(primitive_state, species)
    bulk_velocity = get_bulk_velocity_from_primitive_state(primitive_state)
    momentum_density = bulk_velocity * mass_density
    kinetic_energy_density_val = kinetic_energy_density(mass_density, momentum_density)
    internal_energy_density = get_internal_energy_density_from_primitive_state(primitive_state, species)
    total_energy_density = kinetic_energy_density_val + internal_energy_density
    return construct_conservative_state(mass_density, momentum_density, total_energy_density)
