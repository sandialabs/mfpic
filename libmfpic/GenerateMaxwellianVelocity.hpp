#pragma once

#include <libmfpic/Constants.hpp>

#include <mfem/mfem.hpp>

#include <random>

namespace mfpic {

/**
 * @brief Generate a three-dimensional velocity vector from a drifting Maxwellian distribution.
 *
 * @tparam Generator A UniformRandomBitGenerator type.
 *
 * @param[in]     bulk_velocity              Bulk (mean) velocity of distribution.
 * @param[in]     temperature                Temperature of distribution.
 * @param[in]     mass                       Mass of species.
 * @param[in,out] generator                  A UniformRandomBitGenerator used to generate some random numbers.
 *
 * @returns Velocity randomly drawn the requested drifting Maxwellian distribution.
 */
template <std::uniform_random_bit_generator Generator>
mfem::Vector generateMaxwellianVelocity(
  const mfem::Vector& bulk_velocity,
  double temperature,
  double mass,
  Generator& generator
) {
  assert(bulk_velocity.Size() == 3);
  assert(temperature >= 0.0);
  assert(mass > 0.0);

  const double thermal_speed = sqrt(constants::boltzmann_constant * temperature / mass);
  mfem::Vector velocity = bulk_velocity;
  if (thermal_speed > 0.0) {
    constexpr double peculiar_velocity_distribution_mean = 0.0;
    std::normal_distribution<double> peculiar_velocity_distribution(peculiar_velocity_distribution_mean, thermal_speed);
#ifdef ONE_V
    for (int dimension = 0; dimension < 1; dimension++) {
      velocity[dimension] += peculiar_velocity_distribution(generator);
    }
#else
    for (int dimension = 0; dimension < 3; dimension++) {
      velocity[dimension] += peculiar_velocity_distribution(generator);
    }
#endif
  }

  return velocity;
}

} // namespace mfpic
