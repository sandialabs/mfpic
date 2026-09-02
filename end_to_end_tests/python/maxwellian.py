import numpy as np
from scipy.constants import Boltzmann as boltzmann_constant


def thermal_speed(species, temperature):
    return np.sqrt(boltzmann_constant * temperature / species.mass)


def probability_density_function_1d(velocity, species, bulk_velocity, temperature):
    sigma = thermal_speed(species, temperature)
    inv_sq_sigma = 1.0 / (sigma * sigma)
    velocity_difference = velocity - bulk_velocity
    exponent = inv_sq_sigma * np.power(velocity_difference, 2)
    norm = 1.0 / (np.sqrt(2.0 * np.pi) * sigma)
    probability_density_function = norm * np.exp(-0.5 * exponent)
    return probability_density_function
