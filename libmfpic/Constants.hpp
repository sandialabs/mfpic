#pragma once

#include <cmath>

namespace mfpic {
namespace constants {

/// Exact value [m/s].
static double constexpr speed_of_light = 299792458.0;
/// Exact value [m^2/s^2].
static double constexpr speed_of_light_squared = 89875517873681764.0;
/// Approximate value [H/m]. Error in 10th digit.
static double constexpr permeability = 4 * M_PI * 1.00000000082e-7;
/// Approximate value [F/m]. Same uncertainty as permeability.
static double constexpr permittivity = 1. / (permeability * speed_of_light * speed_of_light);
/// Exact value [J/K].
static double constexpr boltzmann_constant = 1.38064852e-23;
/// Approximate value [kg]. Relative error +/- 12e-9.
static double constexpr atomic_mass_unit = 1.660539040e-27;
/// Exact value [C].
static double constexpr elementary_charge = 1.602176634e-19;
/// Approximate value [kg]. Relative error +/- 5e-10.
static double constexpr electron_mass = 9.10938356e-31;
/// Approximate value [kg]. Relative error in coefficient +/- 9e-11.
static double constexpr proton_mass = 1.007276466879 * atomic_mass_unit;

}
}