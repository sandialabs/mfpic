import numpy as np

speed_of_light = 299792458.0
vacuum_permeability = 4 * np.pi * 1.00000000082e-7
vacuum_permittivity = 1. / (vacuum_permeability * speed_of_light * speed_of_light)
boltzmann_constant = 1.38064852e-23
atomic_mass_unit = 1.660539040e-27
elementary_charge = 1.602176634e-19
electron_mass = 9.10938356e-31
proton_mass = 1.007276466879 * atomic_mass_unit
electron_volt_in_kelvin = elementary_charge / boltzmann_constant