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

  ParticleMoments unrelaxed_particle_moments = particle_operations.getParticleMoments(unrelaxed_particles);
  auto computed_unrelaxed_number_density = unrelaxed_particle_moments.number_density;
  auto computed_unrelaxed_bulk_velocity = unrelaxed_particle_moments.bulk_velocity;
  auto computed_unrelaxed_temperature = unrelaxed_particle_moments.temperature;

  ParticleMoments relaxed_particle_moments = particle_operations.getParticleMoments(relaxed_particles);
  auto computed_relaxed_number_density = relaxed_particle_moments.number_density;
  auto computed_relaxed_bulk_velocity = relaxed_particle_moments.bulk_velocity;
  auto computed_relaxed_temperature = relaxed_particle_moments.temperature;


  EXPECT_EQ(unrelaxed_particles.numParticles(), relaxed_particles.numParticles());
  constexpr double relative_tolerance = 0.1;
  double unrelaxed_number_density = computed_relaxed_number_density.at(default_species)[0];
  double relaxed_number_density = computed_unrelaxed_number_density.at(default_species)[0];
  EXPECT_NEAR(relaxed_number_density, unrelaxed_number_density, relative_tolerance*unrelaxed_number_density);
  for (int idim = 0; idim < 3; idim++) {
    double unrelaxed_bulk_velocity = computed_unrelaxed_bulk_velocity.at(default_species)(idim, 0);
    double relaxed_bulk_velocity = computed_relaxed_bulk_velocity.at(default_species)(idim, 0);
    EXPECT_NEAR(relaxed_bulk_velocity, unrelaxed_bulk_velocity, relative_tolerance*unrelaxed_bulk_velocity);
  }
  double unrelaxed_temperature = computed_unrelaxed_temperature.at(default_species)[0];
  double relaxed_temperature = computed_relaxed_temperature.at(default_species)[0];
  EXPECT_NEAR(relaxed_temperature, unrelaxed_temperature, relative_tolerance*unrelaxed_temperature);
}

double totalKineticEnergy(const ParticleContainer& particles) {
  double total = 0.0;
  for (const Particle& particle : particles) total += 0.5 * particle.weight * particle.species.mass * (particle.velocity * particle.velocity);
  return total;
}

ParticleContainer makeParticlesWithVaryingWeightsInBothElements(int num_particles) {
  ParticleContainer particles;
  for (int iparticle = 0; iparticle < num_particles; iparticle++) {
    const int element = iparticle % num_elems;
    particles.addParticle(Particle{
      .position=mfem::Vector{(element + 0.5) * domain_length / num_elems, 0.0, 0.0},
      .velocity=mfem::Vector{1.0 * iparticle, -2.0 * iparticle + 50.0, 3.0 * (iparticle % 7)},
      .element=element,
      .species=default_species,
      .weight=1.0 + 0.1 * (iparticle % 5)});
  }
  return particles;
}

TEST(BGKRelaxationOperations, LoneRelaxingParticleInElementKeepsItsVelocity) {
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );
  ParticleContainer unrelaxed_particles;
  unrelaxed_particles.addParticle(Particle{.velocity=mfem::Vector{100, 200, 300}, .element=0, .species=default_species, .weight=1.0});
  unrelaxed_particles.addParticle(Particle{.velocity=mfem::Vector{-50, 20, 10}, .element=1, .species=default_species, .weight=1.0});
  unrelaxed_particles.addParticle(Particle{.velocity=mfem::Vector{70, 80, -90}, .element=1, .species=default_species, .weight=2.0});

  ParticleContainer relaxed_particles = unrelaxed_particles;
  constexpr double collision_frequency = 1.0;
  BGKRelaxationOperations relaxer(collision_frequency, default_species);
  constexpr double dt = 1.0e15 / collision_frequency;
  RandomNumberGenerator generator;
  relaxer.performCollisions(dt, generator, relaxed_particles, particle_operations);

  const Particle& unrelaxed_lone_particle = *unrelaxed_particles.begin();
  const Particle& relaxed_lone_particle = *relaxed_particles.begin();
  for (int idim = 0; idim < 3; idim++) {
    EXPECT_THAT(relaxed_lone_particle.velocity[idim], DoubleEq(unrelaxed_lone_particle.velocity[idim]));
  }
}

TEST(BGKRelaxationOperations, EnergyIsConservedOnAverageWithoutExactConservation) {
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );
  const ParticleContainer unrelaxed_particles = makeParticlesWithVaryingWeightsInBothElements(2000);
  const double unrelaxed_energy = totalKineticEnergy(unrelaxed_particles);

  constexpr double collision_frequency = 1.0;
  BGKRelaxationOperations relaxer(collision_frequency, default_species);
  constexpr double dt = 1.0e15 / collision_frequency;
  constexpr int num_trials = 4000;
  double mean_relative_energy_change = 0.0;
  for (int trial = 0; trial < num_trials; trial++) {
    ParticleContainer relaxed_particles = unrelaxed_particles;
    RandomNumberGenerator generator(trial);
    relaxer.performCollisions(dt, generator, relaxed_particles, particle_operations);
    mean_relative_energy_change += (totalKineticEnergy(relaxed_particles) - unrelaxed_energy) / unrelaxed_energy / num_trials;
  }
  EXPECT_NEAR(mean_relative_energy_change, 0.0, 1.0e-2);
}

} // namespace
