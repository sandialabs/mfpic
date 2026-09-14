#include <libmfpic/BGKRelaxationOperations.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/ParticleOperations.hpp>
#include <libmfpic/ReflectingParticleBoundary.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <mfem/mfem.hpp>

#include <ranges>

namespace {

using namespace mfpic;
using namespace ::testing;

const Species default_species{.charge = 1.0, .mass = 1.0};
const std::unordered_map<std::string, Species> one_species {{"one", default_species}};
const std::vector<std::shared_ptr<ParticleBoundaryFactory>> empty_particle_boundary_factory_list;
const std::shared_ptr<ParticleBoundaryFactory> default_reflecting_particle_boundary_factory
  = std::make_shared<ReflectingParticleBoundaryFactory>();
constexpr int num_elems = 2;
constexpr double domain_length = 1.0;
constexpr int order = 1;
mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems, domain_length);

TEST(BGKRelaxationOperations, NothingHappensToParticlesWhenTimestepIsSmallComparedToMeanCollisionTime) {
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );
  ParticleContainer unrelaxed_particles;
  unrelaxed_particles.addParticle(Particle{.velocity=mfem::Vector{100, 200, 300}, .species=default_species, .weight=1.0});
  unrelaxed_particles.addParticle(Particle{.velocity=mfem::Vector{142, 834, 123}, .species=default_species, .weight=1.0});
  unrelaxed_particles.addParticle(Particle{.velocity=mfem::Vector{721, 175, 435}, .species=default_species, .weight=1.0});

  ParticleContainer relaxed_particles = unrelaxed_particles;
  constexpr double collision_frequency = 1.0;
  BGKRelaxationOperations relaxer(collision_frequency, default_species);
  constexpr double dt = 1.0e-15 / collision_frequency;
  RandomNumberGenerator generator;
  relaxer.performCollisions(dt, generator, relaxed_particles, particle_operations);

  ASSERT_EQ(unrelaxed_particles.numParticles(), relaxed_particles.numParticles());
  for (const auto [unrelaxed_particle, relaxed_particle] : std::views::zip(unrelaxed_particles, relaxed_particles)) {
    const mfem::Vector unrelaxed_velocity = unrelaxed_particle.velocity;
    const mfem::Vector relaxed_velocity = relaxed_particle.velocity;
    for (int idim = 0; idim < 3; idim++) {
      EXPECT_THAT(unrelaxed_velocity[idim], DoubleEq(relaxed_velocity[idim]));
    }
  }
}

TEST(BGKRelaxationOperations, EveryParticleRelaxesWhenTimestepIsLargeComparedToMeanCollisionTime) {
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );
  ParticleContainer unrelaxed_particles;
  unrelaxed_particles.addParticle(Particle{.velocity=mfem::Vector{100, 200, 300}, .species=default_species, .weight=1.0});
  unrelaxed_particles.addParticle(Particle{.velocity=mfem::Vector{142, 834, 123}, .species=default_species, .weight=1.0});
  unrelaxed_particles.addParticle(Particle{.velocity=mfem::Vector{721, 175, 435}, .species=default_species, .weight=1.0});

  ParticleContainer relaxed_particles = unrelaxed_particles;
  constexpr double collision_frequency = 1.0;
  BGKRelaxationOperations relaxer(collision_frequency, default_species);
  constexpr double dt = 1.0e15 / collision_frequency;
  RandomNumberGenerator generator;
  relaxer.performCollisions(dt, generator, relaxed_particles, particle_operations);

  ASSERT_EQ(unrelaxed_particles.numParticles(), relaxed_particles.numParticles());
  for (const auto [unrelaxed_particle, relaxed_particle] : std::views::zip(unrelaxed_particles, relaxed_particles)) {
    const mfem::Vector unrelaxed_velocity = unrelaxed_particle.velocity;
    const mfem::Vector relaxed_velocity = relaxed_particle.velocity;
    for (int idim = 0; idim < 3; idim++) {
      EXPECT_THAT(unrelaxed_velocity[idim], Not(DoubleEq(relaxed_velocity[idim])));
    }
  }
}

TEST(BGKRelaxationOperations, FractionOfRelaxedParticlesIsAsExpected) {
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );
  constexpr int num_particles = 1000;
  ParticleContainer unrelaxed_particles;
  for (int iparticle = 0; iparticle < num_particles; iparticle++) {
    unrelaxed_particles.addParticle(Particle{.velocity=mfem::Vector{iparticle, 2*iparticle, 300*iparticle}, .species=default_species, .weight=1.0});
  }

  ParticleContainer relaxed_particles = unrelaxed_particles;
  constexpr double collision_frequency = 1.0;
  BGKRelaxationOperations relaxer(collision_frequency, default_species);
  constexpr double dt = 1.0 / collision_frequency;
  RandomNumberGenerator generator;
  relaxer.performCollisions(dt, generator, relaxed_particles, particle_operations);

  ASSERT_EQ(unrelaxed_particles.numParticles(), relaxed_particles.numParticles());
  int num_relaxed_particles = 0;
  for (const auto [unrelaxed_particle, relaxed_particle] : std::views::zip(unrelaxed_particles, relaxed_particles)) {
    const mfem::Vector unrelaxed_velocity = unrelaxed_particle.velocity;
    const mfem::Vector relaxed_velocity = relaxed_particle.velocity;
    bool particle_is_untouched = true;
    for (int idim = 0; idim < 3; idim++) {
      particle_is_untouched = particle_is_untouched and unrelaxed_velocity[idim] == relaxed_velocity[idim];
    }
    if (not particle_is_untouched) {
      num_relaxed_particles += 1;
    }
  }
  const double relaxation_fraction = -std::expm1(-collision_frequency * dt);
  const double expected_num_relaxed_particles = num_particles * relaxation_fraction;
  const double tolerance = 0.1 * expected_num_relaxed_particles;
  EXPECT_NEAR(num_relaxed_particles, expected_num_relaxed_particles, tolerance);
}

TEST(BGKRelaxationOperations, MomentsAreApproximatelyConservedWhenAllParticlesRelax) {
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );
  constexpr int num_particles = 1000;
  ParticleContainer unrelaxed_particles;
  for (int iparticle = 0; iparticle < num_particles; iparticle++) {
    unrelaxed_particles.addParticle(Particle{.velocity=mfem::Vector{iparticle, 2*iparticle, 3*iparticle}, .species=default_species, .weight=1.0});
  }

  ParticleContainer relaxed_particles = unrelaxed_particles;
  constexpr double collision_frequency = 1.0;
  BGKRelaxationOperations relaxer(collision_frequency, default_species);
  constexpr double dt = 1.0e15 / collision_frequency;
  RandomNumberGenerator generator;
  relaxer.performCollisions(dt, generator, relaxed_particles, particle_operations);

  EXPECT_EQ(unrelaxed_particles.numParticles(), relaxed_particles.numParticles());
  constexpr double relative_tolerance = 0.1;
  double unrelaxed_number_density = particle_operations.getNumberDensity(unrelaxed_particles).at(default_species)[0];
  double relaxed_number_density = particle_operations.getNumberDensity(relaxed_particles).at(default_species)[0];
  EXPECT_NEAR(relaxed_number_density, unrelaxed_number_density, relative_tolerance*unrelaxed_number_density);
  for (int idim = 0; idim < 3; idim++) {
    double unrelaxed_bulk_velocity = particle_operations.getBulkVelocity(unrelaxed_particles).at(default_species)(idim, 0);
    double relaxed_bulk_velocity = particle_operations.getBulkVelocity(relaxed_particles).at(default_species)(idim, 0);
    EXPECT_NEAR(relaxed_bulk_velocity, unrelaxed_bulk_velocity, relative_tolerance*unrelaxed_bulk_velocity);
  }
  double unrelaxed_temperature = particle_operations.getTemperature(unrelaxed_particles).at(default_species)[0];
  double relaxed_temperature = particle_operations.getTemperature(relaxed_particles).at(default_species)[0];
  EXPECT_NEAR(relaxed_temperature, unrelaxed_temperature, relative_tolerance*unrelaxed_temperature);
}

} // namespace
