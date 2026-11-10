#include <libmfpic/AbsorbingParticleBoundary.hpp>
#include <libmfpic/Constants.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/MaxwellianParticleBoundarySource.hpp>
#include <libmfpic/ReflectingParticleBoundary.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

const Species proton_species{.charge = constants::elementary_charge, .mass = constants::proton_mass, .name = "proton"};
const std::unordered_map<std::string, Species> species_map{{"proton", proton_species}};
const std::vector<std::shared_ptr<ParticleBoundaryFactory>> no_particle_boundary_factories;

constexpr double flux_density = 1.0e20;
constexpr double temperature = 11600.0;
constexpr int left_boundary_attribute_1d = 1;

struct BoundarySourceTestSetup {
  std::shared_ptr<mfem::Mesh> mesh;
  Discretization discretization;
  ParticleOperations particle_operations;

  BoundarySourceTestSetup(
    std::shared_ptr<mfem::Mesh> mesh_in,
    std::shared_ptr<ParticleBoundaryFactory> default_particle_boundary_factory)
    : mesh(mesh_in)
    , discretization(mesh.get(), 1)
    , particle_operations(discretization, no_particle_boundary_factories, default_particle_boundary_factory, species_map, 3)
  {}
};

double thermalSpeed() {
  return std::sqrt(constants::boltzmann_constant * temperature / proton_species.mass);
}

TEST(MaxwellianParticleBoundarySource, CreatesRequestedNumberOfParticlesWithFluxWeight) {
  constexpr double domain_length = 1.0;
  constexpr int num_elements = 10;
  BoundarySourceTestSetup setup(
    std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elements, domain_length)),
    std::make_shared<ReflectingParticleBoundaryFactory>());
  constexpr int num_particles = 50;
  MaxwellianParticleBoundarySource source(
    setup.mesh, left_boundary_attribute_1d, proton_species, num_particles, flux_density, temperature, 3);
  RandomNumberGenerator generator;
  const double dt = 1.0e-3 * domain_length / thermalSpeed();

  const ParticleContainer particles = source.loadParticles(generator, setup.particle_operations, dt);

  ASSERT_EQ(particles.numParticles(), num_particles);
  const double expected_weight = flux_density * 1.0 * dt / num_particles;
  for (const Particle& particle : particles) {
    EXPECT_DOUBLE_EQ(particle.weight, expected_weight);
    EXPECT_EQ(particle.species, proton_species);
    EXPECT_TRUE(particle.is_alive);
    EXPECT_GT(particle.particle_distribution_function_value, 0.0);
  }
}

TEST(MaxwellianParticleBoundarySource, ParticlesArePushedInwardByLessThanOneTimestep) {
  constexpr double domain_length = 1.0;
  constexpr int num_elements = 10;
  BoundarySourceTestSetup setup(
    std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elements, domain_length)),
    std::make_shared<ReflectingParticleBoundaryFactory>());
  constexpr int num_particles = 500;
  MaxwellianParticleBoundarySource source(
    setup.mesh, left_boundary_attribute_1d, proton_species, num_particles, flux_density, temperature, 3);
  RandomNumberGenerator generator;
  const double dt = 0.05 * domain_length / thermalSpeed();

  const ParticleContainer particles = source.loadParticles(generator, setup.particle_operations, dt);

  int num_particles_strictly_inside = 0;
  for (const Particle& particle : particles) {
    EXPECT_GT(particle.velocity[0], 0.0);
    EXPECT_GE(particle.position[0], 0.0);
    EXPECT_LE(particle.position[0], particle.velocity[0] * dt);
    const double element_width = domain_length / num_elements;
    EXPECT_EQ(particle.element, std::min(static_cast<int>(particle.position[0] / element_width), num_elements - 1));
    if (particle.position[0] > 0.0) num_particles_strictly_inside++;
  }
  EXPECT_EQ(num_particles_strictly_inside, num_particles);
}

TEST(MaxwellianParticleBoundarySource, ParticlesCrossingAbsorbingBoundaryDuringPartialPushAreRemoved) {
  constexpr double domain_length = 1.0;
  constexpr int num_elements = 4;
  BoundarySourceTestSetup setup(
    std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elements, domain_length)),
    std::make_shared<AbsorbingParticleBoundaryFactory>());
  constexpr int num_particles = 1000;
  MaxwellianParticleBoundarySource source(
    setup.mesh, left_boundary_attribute_1d, proton_species, num_particles, flux_density, temperature, 3);
  RandomNumberGenerator generator;
  const double dt = 2.0 * domain_length / thermalSpeed();

  const ParticleContainer particles = source.loadParticles(generator, setup.particle_operations, dt);

  EXPECT_LT(particles.numParticles(), num_particles);
  EXPECT_GT(particles.numParticles(), 0);
  for (const Particle& particle : particles) {
    EXPECT_TRUE(particle.is_alive);
    EXPECT_LT(particle.position[0], domain_length);
  }
}

TEST(MaxwellianParticleBoundarySource, ParticlesEnter2DDomainThroughTopSide) {
  constexpr int num_elements_per_side = 4;
  std::shared_ptr<mfem::Mesh> mesh = std::make_shared<mfem::Mesh>(
    mfem::Mesh::MakeCartesian2D(num_elements_per_side, num_elements_per_side, mfem::Element::TRIANGLE));
  constexpr int top_boundary_attribute_2d = 3;
  BoundarySourceTestSetup setup(mesh, std::make_shared<ReflectingParticleBoundaryFactory>());
  constexpr int num_particles = 300;
  MaxwellianParticleBoundarySource source(
    setup.mesh, top_boundary_attribute_2d, proton_species, num_particles, flux_density, temperature, 3);
  RandomNumberGenerator generator;
  const double dt = 0.02 / thermalSpeed();

  const ParticleContainer particles = source.loadParticles(generator, setup.particle_operations, dt);

  ASSERT_EQ(particles.numParticles(), num_particles);
  for (const Particle& particle : particles) {
    EXPECT_LT(particle.velocity[1], 0.0);
    EXPECT_LE(particle.position[1], 1.0);
    EXPECT_GE(particle.position[1], 1.0 + particle.velocity[1] * dt * (1.0 + 1.0e-12));
  }
}

TEST(MaxwellianParticleBoundarySource, ThrowsIfBoundaryAttributeHasNoFaces) {
  const auto mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(4, 1.0));
  constexpr int nonexistent_boundary_attribute = 7;

  EXPECT_THROW(
    MaxwellianParticleBoundarySource(mesh, nonexistent_boundary_attribute, proton_species, 1, flux_density, temperature, 3),
    std::logic_error);
}

TEST(MaxwellianParticleBoundarySource, DistributionFunctionIntegratesToNearWallDensityOverInwardVelocities) {
  const auto mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(4, 1.0));
  const MaxwellianParticleBoundarySource source(
    mesh, left_boundary_attribute_1d, proton_species, 1, flux_density, temperature, 1);

  constexpr int num_intervals = 4000;
  const double max_speed = 10.0 * thermalSpeed();
  const double dv = max_speed / num_intervals;
  double integral = 0.0;
  mfem::Vector velocity(1);
  for (int i = 0; i < num_intervals; i++) {
    velocity[0] = (i + 0.5) * dv;
    integral += source.distributionFunctionValue(velocity) * dv;
  }

  const double expected_near_wall_density = flux_density * std::sqrt(M_PI / 2.0) / thermalSpeed();
  constexpr double relative_tolerance = 1.0e-6;
  EXPECT_NEAR(integral, expected_near_wall_density, relative_tolerance * expected_near_wall_density);
}

TEST(MaxwellianParticleBoundarySource, InwardVelocityTimesDistributionFunctionIntegratesToFluxDensityOverInwardVelocities) {
  const auto mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(4, 1.0));
  const MaxwellianParticleBoundarySource source(
    mesh, left_boundary_attribute_1d, proton_species, 1, flux_density, temperature, 1);

  constexpr int num_intervals = 4000;
  const double max_speed = 10.0 * thermalSpeed();
  const double dv = max_speed / num_intervals;
  double integral = 0.0;
  mfem::Vector velocity(1);
  for (int i = 0; i < num_intervals; i++) {
    velocity[0] = (i + 0.5) * dv;
    integral += velocity[0] * source.distributionFunctionValue(velocity) * dv;
  }

  constexpr double relative_tolerance = 1.0e-6;
  EXPECT_NEAR(integral, flux_density, relative_tolerance * flux_density);
}

} // namespace
