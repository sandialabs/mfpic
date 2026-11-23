import numpy as np
from constants import vacuum_permittivity, boltzmann_constant


def plasma_frequency(species, number_density):
    return np.sqrt(number_density * species.charge**2 / (species.mass * vacuum_permittivity))


def debye_length(species, number_density, temperature):
    return np.sqrt(vacuum_permittivity * boltzmann_constant * temperature / (number_density * species.charge**2))