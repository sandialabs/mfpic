#include <libmfpic/GenerateMaxwellianVelocity.hpp>

#include <numbers>

namespace mfpic {

mfem::Vector generateMaxwellianVelocity(
  const mfem::Vector& bulk_velocity,
  double temperature,
  double mass,
  RandomNumberGenerator& generator,
  const int velocity_dims
) {
  assert(bulk_velocity.Size() == 3);
  assert(temperature >= 0.0);
  assert(mass > 0.0);

  const double thermal_speed = sqrt(constants::boltzmann_constant * temperature / mass);
  mfem::Vector velocity = bulk_velocity;
  if (thermal_speed > 0.0) {
    constexpr double peculiar_velocity_distribution_mean = 0.0;
    std::normal_distribution<double> peculiar_velocity_distribution(peculiar_velocity_distribution_mean, thermal_speed);

    for (int dimension = 0; dimension < velocity_dims; dimension++)
      velocity[dimension] += peculiar_velocity_distribution(generator);
  }

  return velocity;
}

mfem::Vector generateMaxwellianFluxVelocity(
  const mfem::Vector& unit_normal,
  double temperature,
  double mass,
  RandomNumberGenerator& generator,
  int velocity_dims
) {
  assert(unit_normal.Size() == 3);
  assert(temperature > 0.0);

  const mfem::Vector zero_bulk_velocity({0.0, 0.0, 0.0});
  mfem::Vector velocity = generateMaxwellianVelocity(zero_bulk_velocity, temperature, mass, generator, velocity_dims);

  const double thermal_speed = sqrt(constants::boltzmann_constant * temperature / mass);
  constexpr double weibull_a_to_get_rayleigh_distribution = 2;
  const double weibull_b_to_get_rayleigh_distribution = std::numbers::sqrt2 * thermal_speed;
  std::weibull_distribution<double> normal_speed_distribution(
    weibull_a_to_get_rayleigh_distribution,
    weibull_b_to_get_rayleigh_distribution);
  const double normal_speed = normal_speed_distribution(generator);

  velocity.Add(normal_speed - velocity * unit_normal, unit_normal);

  return velocity;
}

} // namespace mfpic
