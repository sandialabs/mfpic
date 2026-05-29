#pragma once

#include <libmfpic/Constants.hpp>

#include <mfem/mfem.hpp>

#include <random>

namespace mfpic {

/**
 * @brief Generate a three-dimensional velocity vector from the product of 1D Kappa distributions.
 *
 * @tparam Generator A UniformRandomBitGenerator type.
 *
 * @param[in]     bulk_velocity              Bulk (mean) velocity of distribution.
 * @param[in]     temperature                Temperature of distribution.
 * @param[in]     kappa                      Kappa parameter in the kappa-distribution.
 * @param[in]     mass                       Mass of species.
 * @param[in,out] generator                  A UniformRandomBitGenerator used to generate some random numbers.
 *
 * @returns Velocity randomly drawn from the requested Kappa distribution.
 */
template <std::uniform_random_bit_generator Generator>
mfem::Vector generateKappaVelocity(
  const mfem::Vector& bulk_velocity,
  double temperature,
  double kappa,
  double mass,
  Generator& generator
) {
  assert(bulk_velocity.Size() == 3);
  assert(temperature >= 0.0);
  assert(mass > 0.0);

  const double v_thermal_squared = constants::boltzmann_constant * temperature / mass;
  mfem::Vector velocity = bulk_velocity;
  if (v_thermal_squared > 0.0) {
    const double nu = 2.0 * kappa + 1.0;
    const double scale = std::sqrt(kappa * 2.0 * v_thermal_squared / nu);
    std::student_t_distribution<double> velocity_distribution(nu);
    for (int dimension = 0; dimension < 3; dimension++) {
      velocity[dimension] += scale * velocity_distribution(generator);
    }
  }

  return velocity;
}

/**
 * @brief Generate a three-dimensional velocity vector from an isotropic Kappa distribution.
 *
 * @tparam Generator A UniformRandomBitGenerator type.
 *
 * @param[in]     bulk_velocity              Bulk (mean) velocity of distribution.
 * @param[in]     temperature                Temperature of distribution.
 * @param[in]     kappa                      Kappa parameter in the kappa-distribution.
 * @param[in]     mass                       Mass of species.
 * @param[in,out] generator                  A UniformRandomBitGenerator used to generate some random numbers.
 *
 * @returns Velocity randomly drawn from the requested Kappa distribution.
 */
template <std::uniform_random_bit_generator Generator>
mfem::Vector generateIsotropicKappaVelocity(
  const mfem::Vector& bulk_velocity,
  double temperature,
  double kappa,
  double mass,
  Generator& generator)
{
  assert(bulk_velocity.Size() == 3);
  assert(temperature >= 0.0);
  assert(mass > 0.0);
  assert(kappa > 1.5);

  mfem::Vector velocity = bulk_velocity;

  if (temperature == 0.0) {
    return velocity;
  }

  const double nu = 2.0 * kappa - 1.0;
  const double w_squared =
  ((nu - 2.0) / nu) *
  constants::boltzmann_constant *
  temperature / mass;
  const double w = std::sqrt(w_squared);
  std::normal_distribution<double> normal_distribution(0.0, 1.0);
  std::gamma_distribution<double> gamma_distribution(
    0.5 * nu,
    2.0);
  const double chi_squared = gamma_distribution(generator);
  const double scale = std::sqrt(nu / chi_squared);
  for (int dimension = 0; dimension < 3; ++dimension) {
    velocity(dimension) +=
      w * scale * normal_distribution(generator);
  }
  return velocity;
}

} // namespace mfpic
