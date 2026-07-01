#include <libmfpic/GenerateKappaVelocity.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

constexpr double mass = 1.0;

TEST(GenerateKappaVelocity, GeneratedVelocityEqualsBulkVelocityWhenTemperatureIsZero) {
  const mfem::Vector bulk_velocity({1.0, 2.0, 3.0});
  constexpr double temperature = 0.0;
  constexpr double kappa = 2.0;
  RandomNumberGenerator generator;

  mfem::Vector velocity = generateKappaVelocity(bulk_velocity, temperature, kappa, mass, generator);

  for (int dimension = 0; dimension < 3; dimension++) {
    EXPECT_DOUBLE_EQ(velocity[dimension], bulk_velocity[dimension]);
  }
}

TEST(GenerateKappaVelocity, GeneratedVelocitiesAreMaxwellianWithLargeKappa) {
  const mfem::Vector nominal_bulk_velocity({300.0, 600.0, 1000.0});
  constexpr double temperature = 11600.0;
  constexpr int num_samples = 20000;
  constexpr double kappa = 1e16;
  RandomNumberGenerator generator;

  std::vector<mfem::Vector> generated_velocities;
  for (int i = 0; i < num_samples; i++) {
    generated_velocities.emplace_back(generateKappaVelocity(nominal_bulk_velocity, temperature, kappa, mass, generator));
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

TEST(GenerateKappaVelocity, VelocitiesMeanAndStdAreCorrect) {
  const mfem::Vector nominal_bulk_velocity({300.0, 600.0, 1000.0});
  constexpr double temperature = 11600.0;
  constexpr int num_samples = 20000;
  constexpr double kappa = 2.0;
  RandomNumberGenerator generator;

  std::vector<mfem::Vector> generated_velocities;
  for (int i = 0; i < num_samples; i++) {
    generated_velocities.emplace_back(generateKappaVelocity(nominal_bulk_velocity, temperature, kappa, mass, generator));
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
  constexpr double expected_sample_variance_maxwellian = 3.0 * constants::boltzmann_constant * temperature / mass;
  constexpr double expected_sample_variance_kappa = expected_sample_variance_maxwellian * (kappa / (kappa - 0.5));
  EXPECT_NEAR(sample_variance, expected_sample_variance_kappa, relative_tolerance * expected_sample_variance_kappa);
}

TEST(GenerateIsotropicKappaVelocity, VelocitiesMeanAndVarianceAreCorrect)
{
  const mfem::Vector nominal_bulk_velocity(
    {300.0, 600.0, 1000.0});

  constexpr double temperature = 11600.0;
  constexpr double kappa = 2.0;
  constexpr int num_samples = 20000;

  RandomNumberGenerator generator;

  std::vector<mfem::Vector> generated_velocities;

  for (int i = 0; i < num_samples; ++i) {
    generated_velocities.emplace_back(
      generateIsotropicKappaVelocity(
        nominal_bulk_velocity,
        temperature,
        kappa,
        mass,
        generator));
  }

  // Compute sample mean
  mfem::Vector actual_bulk_velocity({0.0, 0.0, 0.0});

  for (const auto& velocity : generated_velocities) {
    actual_bulk_velocity += velocity;
  }

  actual_bulk_velocity /= num_samples;

  constexpr double relative_tolerance = 0.1;

  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(
      actual_bulk_velocity[i],
      nominal_bulk_velocity[i],
      relative_tolerance *
      nominal_bulk_velocity[i]);
  }

  double sample_variance = 0.0;

  for (const auto& velocity : generated_velocities) {

    mfem::Vector relative_velocity = velocity;
    relative_velocity -= actual_bulk_velocity;

    sample_variance +=
      relative_velocity * relative_velocity;
  }

  sample_variance /= (num_samples - 1);

  const double expected_sample_variance =
    3.0 *
    constants::boltzmann_constant *
    temperature / mass;

  EXPECT_NEAR(
    sample_variance,
    expected_sample_variance,
    relative_tolerance *
    expected_sample_variance);
}

} // namespace
