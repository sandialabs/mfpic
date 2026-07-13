#pragma once

#include <libmfpic/Constants.hpp>
#include <libmfpic/RandomNumberGenerator.hpp>

#include <mfem/mfem.hpp>

namespace mfpic {

/**
 * @brief Generate a three-dimensional velocity vector from a drifting Maxwellian distribution.
 *
 * @param[in]     bulk_velocity              Bulk (mean) velocity of distribution.
 * @param[in]     temperature                Temperature of distribution.
 * @param[in]     mass                       Mass of species.
 * @param[in,out] generator                  Random number generator.
 *
 * @returns Velocity randomly drawn the requested drifting Maxwellian distribution.
 */
mfem::Vector generateMaxwellianVelocity(
  const mfem::Vector& bulk_velocity,
  double temperature,
  double mass,
  RandomNumberGenerator& generator
);

mfem::Vector generate1DMaxwellianVelocity(
  const mfem::Vector& bulk_velocity,
  double temperature,
  double mass,
  RandomNumberGenerator& generator
);

} // namespace mfpic
