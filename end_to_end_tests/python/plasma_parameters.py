import numpy as np
from scipy.constants import epsilon_0 as vacuum_permittivity
from scipy.constants import Boltzmann as boltzmann_constant


def plasma_frequency(species, number_density):
    return np.sqrt(number_density * species.charge * species.charge_over_mass / vacuum_permittivity)


def plasma_period(species, number_density):
    return 2.0 * np.pi / plasma_frequency(species, number_density)


def debye_length(species, number_density, temperature):
    return np.sqrt(
        vacuum_permittivity * boltzmann_constant * temperature /
        (number_density * np.power(species.charge, 2))
    )