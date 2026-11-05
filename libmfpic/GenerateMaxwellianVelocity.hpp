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
  RandomNumberGenerator& generator,
  const int velocity_dims=3
);

/**
 * @brief Generate the velocity of a particle crossing a surface from a non-drifting Maxwellian gas.
 *
 * The component along @a unit_normal is drawn from the flux-weighted (Rayleigh) distribution
 * \f$ f(v_n) \propto v_n \exp(-v_n^2 / 2\sigma^2) \f$ for \f$ v_n = \vec{v} \cdot \hat{n} > 0 \f$,
 * where \f$ \sigma = \sqrt{k_B T / m} \f$.
 * The components perpendicular to @a unit_normal are drawn from a Maxwellian.
 *
 * @param[in]     unit_normal        Unit normal of surface that the particle crosses, of length 3.
 * @param[in]     temperature        Temperature of distribution, must be positive.
 * @param[in]     mass               Mass of species.
 * @param[in,out] generator          Random number generator.
 * @param[in]     velocity_dims      Number of velocity dimensions.
 *
 * @returns Velocity having a positive component along @a unit_normal.
 */
mfem::Vector generateMaxwellianFluxVelocity(
  const mfem::Vector& unit_normal,
  double temperature,
  double mass,
  RandomNumberGenerator& generator,
  int velocity_dims
);

} // namespace mfpic
