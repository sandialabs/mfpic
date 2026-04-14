#include <libmfpic/LoadParticles.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

constexpr Species default_species;
const mfem::Vector zero_vector({0.0, 0.0, 0.0});
const int num_elems = 10;
std::shared_ptr<mfem::Mesh> simple_mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elems));

TEST(LoadParticles, NoParticlesAddedWhenNoParticlesRequested) {
  const SourceStateParameters source_state_parameters{
    .number_density = 1.0e18,
    .bulk_velocity = zero_vector,
    .temperature = 121325.0,
  };
  constexpr int num_particles = 0;
  std::default_random_engine generator;

  ParticleContainer particles = loadParticles(
    ConstantSourceParameters(default_species, source_state_parameters, num_particles),
    generator,
    simple_mesh
  );

  ASSERT_EQ(particles.numParticles(), 0);
}

TEST(LoadParticles, NumLoadedParticlesIsAsRequested) {
  const SourceStateParameters source_state_parameters{
    .number_density = 1.0e18,
    .bulk_velocity = zero_vector,
    .temperature = 12415.0,
  };
  constexpr int num_particles = 40;
  std::default_random_engine generator;

  ParticleContainer particles = loadParticles(
    ConstantSourceParameters(default_species, source_state_parameters, num_particles),
    generator,
    simple_mesh
  );

  ASSERT_EQ(particles.numParticles(), num_particles);
}

TEST(LoadParticles, LoadedParticlesAllUseBulkVelocityWithZeroTemperature) {
  const mfem::Vector bulk_velocity({1.0, 2.0, 3.0});
  const SourceStateParameters source_state_parameters{
    .number_density = 1.0e18,
    .bulk_velocity = bulk_velocity,
    .temperature = 0.0,
    .kappa = 2.0,
  };
  constexpr int num_particles = 5;
  std::default_random_engine generator;

  ParticleContainer particles = loadParticles(
    ConstantSourceParameters(default_species, source_state_parameters, num_particles),
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

TEST(LoadParticles, ParticleWeightSetCorrectly) {
  constexpr double number_density = 1.0e18;
  const SourceStateParameters source_state_parameters{
    .number_density = number_density,
    .bulk_velocity = zero_vector,
    .temperature = 0.0,
  };
  constexpr int num_particles = 20;
  std::default_random_engine generator;

  ParticleContainer particles = loadParticles(
    ConstantSourceParameters(default_species, source_state_parameters, num_particles),
    generator,
    simple_mesh
  );

  constexpr double mesh_volume = 1.0;
  constexpr double expected_weight = number_density * mesh_volume / num_particles;
  for (const Particle& particle : particles) {
    EXPECT_DOUBLE_EQ(particle.weight, expected_weight);
  }
}

TEST(LoadParticles, ParticlesAreUniformlyDistributedInSpace) {
  const SourceStateParameters source_state_parameters{
    .number_density = 1.0e18,
    .bulk_velocity = zero_vector,
    .temperature = 0.0,
  };
  constexpr int num_particles = 20000;
  std::default_random_engine generator;

  ParticleContainer particles = loadParticles(
    ConstantSourceParameters(default_species, source_state_parameters, num_particles),
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

TEST(LoadParticles, ParticleVelocitiesAreMaxwellianWhenMaxwellianParticlesAreRequested) {
  const mfem::Vector nominal_bulk_velocity({300.0, 600.0, 1000.0});
  constexpr double temperature = 11600.0;
  const SourceStateParameters source_state_parameters{
    .number_density = 1.0e18,
    .bulk_velocity = nominal_bulk_velocity,
    .temperature = temperature,
  };
  constexpr int num_particles = 20000;
  std::default_random_engine generator;

  ParticleContainer particles = loadParticles(
    ConstantSourceParameters(default_species, source_state_parameters, num_particles),
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

TEST(LoadParticles, ParticleVelocitiesAreMaxwellianWithLargeKappa) {
  const mfem::Vector nominal_bulk_velocity({300.0, 600.0, 1000.0});
  constexpr double temperature = 11600.0;
  const SourceStateParameters source_state_parameters{
    .number_density = 1.0e18,
    .bulk_velocity = nominal_bulk_velocity,
    .temperature = temperature,
    .kappa = 1e16,
  };
  constexpr int num_particles = 20000;
  std::default_random_engine generator;

  ParticleContainer particles = loadParticles(
    ConstantSourceParameters(default_species, source_state_parameters, num_particles),
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

TEST(LoadParticles, ParticleVelocitiesMeanAndStdAreCorrectWhenKappaDistributionIsRequested) {
  const mfem::Vector nominal_bulk_velocity({300.0, 600.0, 1000.0});
  constexpr double temperature = 11600.0;
  constexpr int num_particles = 20000;
  constexpr double kappa = 2.0;
  const SourceStateParameters source_state_parameters{
    .number_density = 1.0e18,
    .bulk_velocity = nominal_bulk_velocity,
    .temperature = temperature,
    .kappa = kappa,
  };
  std::default_random_engine generator;

  ParticleContainer particles = loadParticles(
    ConstantSourceParameters(default_species, source_state_parameters, num_particles),
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
  constexpr double expected_sample_variance_maxwellian = 3.0 * constants::boltzmann_constant * temperature / default_species.mass;
  constexpr double expected_sample_variance_kappa = expected_sample_variance_maxwellian * (kappa / (kappa - 0.5));
  EXPECT_NEAR(sample_variance, expected_sample_variance_kappa, relative_tolerance * expected_sample_variance_kappa);
}

TEST(LoadUniformMaxwellianParticles, ParticleDistributionFunctionsAreMaxwellian) {
  Species species{.charge = -constants::elementary_charge, .mass = constants::electron_mass};
  constexpr double number_density = 1e22;
  constexpr double temperature = 300;
  mfem::Vector bulk_velocity({1.0,2.0,3.0});
  constexpr int num_particles = 1;
  std::mt19937 generator;

  const SourceStateParameters source_state_parameters{
    .number_density = number_density,
    .bulk_velocity = bulk_velocity,
    .temperature = temperature,
  };

  ParticleContainer particles = loadParticles(
    ConstantSourceParameters(species, source_state_parameters, num_particles),
    generator,
    simple_mesh
  );

  mfem::Vector prim(5);
  prim(euler::PrimitiveVariables::NUMBER_DENSITY) = number_density;
  prim(euler::PrimitiveVariables::X_BULK_VELOCITY) = bulk_velocity(0);
  prim(euler::PrimitiveVariables::Y_BULK_VELOCITY) = bulk_velocity(1);
  prim(euler::PrimitiveVariables::Z_BULK_VELOCITY) = bulk_velocity(2);
  prim(euler::PrimitiveVariables::TEMPERATURE) = temperature;
  for (const Particle& particle : particles) {
    const double expected_particle_distribution_value = euler::evaluateMaxwellian(prim, particle.velocity, species);
    EXPECT_DOUBLE_EQ(particle.particle_distribution_function_value, expected_particle_distribution_value);
  }
}

TEST(LoadParticles, LoadedParticlesRespectSodDiscontinuity) {
  constexpr int num_particles = 100;
  constexpr double left_number_density = 1.0;
  constexpr double discontinuity_location = 0.7;
  const SourceStateParameters left_state_parameters{.number_density = left_number_density};
  const SourceStateParameters right_state_parameters{.number_density = 0.0};
  const SodSourceParameters sod_parameters(
    default_species,
    discontinuity_location,
    left_state_parameters,
    right_state_parameters,
    num_particles
  );
  std::default_random_engine generator;

  ParticleContainer particles = loadParticles(
    sod_parameters,
    generator,
    simple_mesh
  );

  constexpr double expected_physical_particles = discontinuity_location * left_number_density;
  constexpr double expected_particle_weight = expected_physical_particles / num_particles;
  for (const Particle& particle : particles) {
    EXPECT_DOUBLE_EQ(particle.weight, expected_particle_weight);
    EXPECT_LE(particle.position[0], discontinuity_location);
  }
}

} // namespace
