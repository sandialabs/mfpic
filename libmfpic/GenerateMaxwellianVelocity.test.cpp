#include <libmfpic/GenerateMaxwellianVelocity.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

using namespace mfpic;

constexpr double mass = 1.0;

TEST(GenerateMaxwellianVelocity, GeneratedVelocityEqualsBulkVelocityWhenTemperatureIsZero) {
  const mfem::Vector bulk_velocity({1.0, 2.0, 3.0});
  constexpr double temperature = 0.0;
  RandomNumberGenerator generator;

  mfem::Vector velocity = generateMaxwellianVelocity(bulk_velocity, temperature, mass, generator);

  for (int dimension = 0; dimension < 3; dimension++) {
    EXPECT_DOUBLE_EQ(velocity[dimension], bulk_velocity[dimension]);
  }
}

TEST(GenerateMaxwellianVelocity, GeneratedVelocitiesAreMaxwellian) {
  const mfem::Vector nominal_bulk_velocity({300.0, 600.0, 1000.0});
  constexpr double temperature = 11600.0;
  constexpr int num_samples = 20000;
  RandomNumberGenerator generator;

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

struct FluxVelocityStatistics {
  double min_normal_speed;
  double mean_normal_speed;
  double mean_squared_normal_speed;
  mfem::Vector mean_tangential_velocity;
  double mean_squared_tangential_speed;
};

FluxVelocityStatistics sampleFluxVelocityStatistics(
  const mfem::Vector& unit_normal,
  double temperature,
  int num_samples,
  int velocity_dims
) {
  RandomNumberGenerator generator;
  FluxVelocityStatistics statistics{
    .min_normal_speed = std::numeric_limits<double>::max(),
    .mean_normal_speed = 0.0,
    .mean_squared_normal_speed = 0.0,
    .mean_tangential_velocity = mfem::Vector({0.0, 0.0, 0.0}),
    .mean_squared_tangential_speed = 0.0,
  };
  for (int i = 0; i < num_samples; i++) {
    const mfem::Vector velocity = generateMaxwellianFluxVelocity(
      unit_normal,
      temperature,
      mass,
      generator,
      velocity_dims
    );
    const double normal_speed = velocity * unit_normal;
    mfem::Vector tangential_velocity = velocity;
    tangential_velocity.Add(-normal_speed, unit_normal);

    statistics.min_normal_speed = std::min(statistics.min_normal_speed, normal_speed);
    statistics.mean_normal_speed += normal_speed / num_samples;
    statistics.mean_squared_normal_speed += normal_speed * normal_speed / num_samples;
    statistics.mean_tangential_velocity.Add(1.0 / num_samples, tangential_velocity);
    statistics.mean_squared_tangential_speed += (tangential_velocity * tangential_velocity) / num_samples;
  }
  return statistics;
}

void testFluxVelocitiesAreFluxWeightedHalfMaxwellian(const mfem::Vector& inward_unit_normal, int velocity_dims) {
  constexpr double temperature = 11600.0;
  constexpr int num_samples = 40000;
  constexpr double relative_tolerance = 0.05;
  const double thermal_speed = std::sqrt(constants::boltzmann_constant * temperature / mass);
  const int num_tangential_dims = velocity_dims - 1;

  const FluxVelocityStatistics statistics = sampleFluxVelocityStatistics(
    inward_unit_normal, temperature, num_samples, velocity_dims);

  EXPECT_GT(statistics.min_normal_speed, 0.0);

  const double expected_mean_normal_speed = thermal_speed * std::sqrt(M_PI / 2.0);
  EXPECT_NEAR(statistics.mean_normal_speed, expected_mean_normal_speed, relative_tolerance * expected_mean_normal_speed);

  const double expected_mean_squared_normal_speed = 2.0 * thermal_speed * thermal_speed;
  EXPECT_NEAR(
    statistics.mean_squared_normal_speed,
    expected_mean_squared_normal_speed,
    relative_tolerance * expected_mean_squared_normal_speed);

  for (int dimension = 0; dimension < 3; dimension++) {
    EXPECT_NEAR(statistics.mean_tangential_velocity[dimension], 0.0, relative_tolerance * thermal_speed);
  }

  const double expected_mean_squared_tangential_speed = num_tangential_dims * thermal_speed * thermal_speed;
  EXPECT_NEAR(
    statistics.mean_squared_tangential_speed,
    expected_mean_squared_tangential_speed,
    relative_tolerance * std::max(expected_mean_squared_tangential_speed, thermal_speed * thermal_speed));
}

TEST(GenerateMaxwellianVelocity, FluxVelocitiesAreFluxWeightedHalfMaxwellianForAxisAlignedNormal) {
  constexpr int num_velocity_dims = 3;
  testFluxVelocitiesAreFluxWeightedHalfMaxwellian(mfem::Vector({1.0, 0.0, 0.0}), num_velocity_dims);
  testFluxVelocitiesAreFluxWeightedHalfMaxwellian(mfem::Vector({0.0, -1.0, 0.0}), num_velocity_dims);
}

TEST(GenerateMaxwellianVelocity, FluxVelocitiesAreFluxWeightedHalfMaxwellianForObliqueNormal) {
  mfem::Vector oblique_normal({1.0, -2.0, 0.5});
  oblique_normal /= oblique_normal.Norml2();
  constexpr int num_velocity_dims = 3;
  testFluxVelocitiesAreFluxWeightedHalfMaxwellian(oblique_normal, num_velocity_dims);
}

TEST(GenerateMaxwellianVelocity, FluxVelocitiesRespectReducedVelocityDimensions) {
  constexpr int num_velocity_dims = 1;
  testFluxVelocitiesAreFluxWeightedHalfMaxwellian(mfem::Vector({-1.0, 0.0, 0.0}), num_velocity_dims);

  constexpr double temperature = 300.0;
  RandomNumberGenerator generator;
  const mfem::Vector velocity = generateMaxwellianFluxVelocity(
    mfem::Vector({1.0, 0.0, 0.0}),
    temperature,
    mass,
    generator,
    num_velocity_dims);
  EXPECT_DOUBLE_EQ(velocity[1], 0.0);
  EXPECT_DOUBLE_EQ(velocity[2], 0.0);
}

} // namespace
