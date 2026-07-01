#include <libmfpic/GenerateMaxwellianVelocity.hpp>

namespace mfpic {

mfem::Vector generateMaxwellianVelocity(
  const mfem::Vector& bulk_velocity,
  double temperature,
  double mass,
  RandomNumberGenerator& generator
) {
  assert(bulk_velocity.Size() == 3);
  assert(temperature >= 0.0);
  assert(mass > 0.0);

  const double thermal_speed = sqrt(constants::boltzmann_constant * temperature / mass);
  mfem::Vector velocity = bulk_velocity;
  if (thermal_speed > 0.0) {
    constexpr double peculiar_velocity_distribution_mean = 0.0;
    std::normal_distribution<double> peculiar_velocity_distribution(peculiar_velocity_distribution_mean, thermal_speed);
    for (int dimension = 0; dimension < 3; dimension++) {
      velocity[dimension] += peculiar_velocity_distribution(generator);
    }
  }

  return velocity;
}

} // namespace mfpic
