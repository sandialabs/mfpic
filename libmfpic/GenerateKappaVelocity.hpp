#pragma once

#include <libmfpic/Constants.hpp>
#include <libmfpic/RandomNumberGenerator.hpp>

#include <mfem/mfem.hpp>

namespace mfpic {

/**
 * @brief Generate a three-dimensional velocity vector from the product of 1D Kappa distributions.
 *
 * @param[in]     bulk_velocity              Bulk (mean) velocity of distribution.
 * @param[in]     temperature                Temperature of distribution.
 * @param[in]     kappa                      Kappa parameter in the kappa-distribution.
 * @param[in]     mass                       Mass of species.
 * @param[in,out] generator                  Random number generator.
 *
 * @returns Velocity randomly drawn from the requested Kappa distribution.
 */
mfem::Vector generateKappaVelocity(
  const mfem::Vector& bulk_velocity,
  double temperature,
  double kappa,
  double mass,
  RandomNumberGenerator& generator
);

/**
 * @brief Generate a three-dimensional velocity vector from an isotropic Kappa distribution.
 *
 * @param[in]     bulk_velocity              Bulk (mean) velocity of distribution.
 * @param[in]     temperature                Temperature of distribution.
 * @param[in]     kappa                      Kappa parameter in the kappa-distribution.
 * @param[in]     mass                       Mass of species.
 * @param[in,out] generator                  Random number generator.
 *
 * @returns Velocity randomly drawn from the requested Kappa distribution.
 */
mfem::Vector generateIsotropicKappaVelocity(
  const mfem::Vector& bulk_velocity,
  double temperature,
  double kappa,
  double mass,
  RandomNumberGenerator& generator);

} // namespace mfpic
