#include <libmfpic/GenerateMaxwellianVelocity.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

constexpr double mass = 1.0;
const mfem::Vector zero_vector({0.0, 0.0, 0.0});

TEST(GenerateMaxwellianVelocity, GeneratedVelocityEqualsBulkVelocityWhenTemperatureIsZero) {
  const mfem::Vector bulk_velocity({1.0, 2.0, 3.0});
  constexpr double temperature = 0.0;
  std::default_random_engine generator;

  mfem::Vector velocity = generateMaxwellianVelocity(bulk_velocity, temperature, mass, generator);

  for (int dimension = 0; dimension < 3; dimension++) {
    EXPECT_DOUBLE_EQ(velocity[dimension], bulk_velocity[dimension]);
  }
}

TEST(GenerateMaxwellianVelocity, GeneratedVelocitiesAreMaxwellian) {
  const mfem::Vector nominal_bulk_velocity({300.0, 600.0, 1000.0});
  constexpr double temperature = 11600.0;
  constexpr int num_samples = 20000;
  std::default_random_engine generator;

  std::vector<mfem::Vector> generated_velocities;
  for (int i = 0; i < num_samples; i++) {
    generated_velocities.emplace_back(generateMaxwellianVelocity(nominal_bulk_velocity, temperature, mass, generator));
  }

  mfem::Vector actual_bulk_velocity({0.0, 0.0, 0.0});
  for (const mfem::Vector& generated_velocity : generated_velocities) {
    actual_bulk_velocity += generated_velocity;
  }
  actual_bulk_velocity /= num_samples;
  constexpr double relative_tolerance = 0.1;
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(actual_bulk_velocity[i], nominal_bulk_velocity[i], relative_tolerance * nominal_bulk_velocity[i]);
  }

  double sample_variance = 0.0;
  for (const mfem::Vector& generated_velocity : generated_velocities) {
    mfem::Vector relative_velocity = generated_velocity;
    relative_velocity -= actual_bulk_velocity;
    sample_variance += relative_velocity * relative_velocity;
  }
  sample_variance /= num_samples - 1;
  constexpr double expected_sample_variance = 3.0 * constants::boltzmann_constant * temperature / mass;
  EXPECT_NEAR(sample_variance, expected_sample_variance, relative_tolerance * expected_sample_variance);
}

} // namespace
