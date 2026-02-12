#include <libmfpic/LoadUniformMaxwellianParticles.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

constexpr Species default_species;
const mfem::Vector zero_vector({0.0, 0.0, 0.0});
const int num_elems = 10;
std::shared_ptr<mfem::Mesh> simple_mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elems));

TEST(LoadUniformMaxwellianParticles, NoParticlesAddedWhenNoParticlesRequested) {
  constexpr double temperature = 121325.0;
  constexpr double number_density = 1.0e18;
  constexpr int num_particles = 0;
  std::mt19937 generator;

  ParticleContainer particles = loadUniformMaxwellianParticles(
    default_species,
    zero_vector,
    temperature,
    number_density,
    num_particles,
    generator,
    simple_mesh
  );

  ASSERT_EQ(particles.numParticles(), 0);
}

TEST(LoadUniformMaxwellianParticles, NumLoadedParticlesIsAsRequested) {
  constexpr double temperature = 12415.0;
  constexpr double number_density = 1.0e18;
  constexpr int num_particles = 40;
  std::mt19937 generator;

  ParticleContainer particles = loadUniformMaxwellianParticles(
    default_species,
    zero_vector,
    temperature,
    number_density,
    num_particles,
    generator,
    simple_mesh
  );

  ASSERT_EQ(particles.numParticles(), num_particles);
}

TEST(LoadUniformMaxwellianParticles, LoadedParticlesAllUseBulkVelocityWithZeroTemperature) {
  const mfem::Vector bulk_velocity({1.0, 2.0, 3.0});
  constexpr double temperature = 0.0;
  constexpr double number_density = 1.0e18;
  constexpr int num_particles = 5;
  std::mt19937 generator;

  ParticleContainer particles = loadUniformMaxwellianParticles(
    default_species,
    bulk_velocity,
    temperature,
    number_density,
    num_particles,
    generator,
    simple_mesh
  );

  ASSERT_EQ(particles.numParticles(), num_particles);
  for (const Particle& particle : particles) {
    for (int dimension = 0; dimension < 3; dimension++) {
      EXPECT_DOUBLE_EQ(particle.velocity[dimension], bulk_velocity[dimension]);
    }
  }
}

TEST(LoadUniformMaxwellianParticles, ParticleWeightSetCorrectly) {
  constexpr double temperature = 0.0;
  constexpr double number_density = 1.0e18;
  constexpr int num_particles = 20;
  std::mt19937 generator;

  ParticleContainer particles = loadUniformMaxwellianParticles(
    default_species,
    zero_vector,
    temperature,
    number_density,
    num_particles,
    generator,
    simple_mesh
  );

  constexpr double mesh_volume = 1.0;
  constexpr double expected_weight = number_density * mesh_volume / num_particles;
  for (const Particle& particle : particles) {
    EXPECT_DOUBLE_EQ(particle.weight, expected_weight);
  }
}

TEST(LoadUniformMaxwellianParticles, ParticlesAreUniformlyDistributedInSpace) {
  constexpr double temperature = 0.0;
  constexpr double number_density = 1.0e18;
  constexpr int num_particles = 20000;
  std::mt19937 generator;

  ParticleContainer particles = loadUniformMaxwellianParticles(
    default_species,
    zero_vector,
    temperature,
    number_density,
    num_particles,
    generator,
    simple_mesh
  );

  constexpr int num_bins = 5;
  const int num_elements = simple_mesh->GetNE();
  const double dx = 1.0 / num_elements;
  std::vector<int> bins(num_bins, 0);
  for (const Particle& particle : particles) {
    ASSERT_DOUBLE_EQ(particle.position[1], 0.0);
    ASSERT_DOUBLE_EQ(particle.position[2], 0.0);
    const int expected_element = particle.position[0] / dx;
    ASSERT_EQ(particle.element, expected_element);

    bins[particle.position[0] * num_bins] += 1;
  }
  constexpr int expected_particles_generated_per_bin = num_particles / num_bins;
  const double absolute_tolerance = 0.1 * expected_particles_generated_per_bin;
  for (int num_particles_generated_in_bin : bins) {
    EXPECT_NEAR(num_particles_generated_in_bin, expected_particles_generated_per_bin, absolute_tolerance);
  }
}

TEST(LoadUniformMaxwellianParticles, ParticleVelocitiesAreMaxwellian) {
  const mfem::Vector nominal_bulk_velocity({300.0, 600.0, 1000.0});
  constexpr double temperature = 11600.0;
  constexpr double number_density = 1.0e18;
  constexpr int num_particles = 20000;
  std::mt19937 generator;

  ParticleContainer particles = loadUniformMaxwellianParticles(
    default_species,
    nominal_bulk_velocity,
    temperature,
    number_density,
    num_particles,
    generator,
    simple_mesh
  );

  mfem::Vector actual_bulk_velocity({0.0, 0.0, 0.0});
  for (const Particle& particle : particles) {
    actual_bulk_velocity += particle.velocity;
  }
  actual_bulk_velocity /= num_particles;
  constexpr double relative_tolerance = 0.1;
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(actual_bulk_velocity[i], nominal_bulk_velocity[i], relative_tolerance * nominal_bulk_velocity[i]);
  }

  double sample_variance = 0.0;
  for (const Particle& particle : particles) {
    mfem::Vector relative_velocity = particle.velocity;
    relative_velocity -= actual_bulk_velocity;
    sample_variance += relative_velocity * relative_velocity;
  }
  sample_variance /= num_particles - 1;
  constexpr double expected_sample_variance = 3.0 * constants::boltzmann_constant * temperature / default_species.mass;
  EXPECT_NEAR(sample_variance, expected_sample_variance, relative_tolerance * expected_sample_variance);
}

} // namespace
