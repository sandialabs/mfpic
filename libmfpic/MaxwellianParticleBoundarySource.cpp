#include <libmfpic/Constants.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/GenerateMaxwellianVelocity.hpp>
#include <libmfpic/MaxwellianParticleBoundarySource.hpp>

namespace mfpic {

MaxwellianParticleBoundarySource::MaxwellianParticleBoundarySource(
  std::shared_ptr<mfem::Mesh> mesh,
  int boundary_attribute,
  const Species& species,
  int num_particles_per_timestep,
  double flux_density,
  double temperature,
  int velocity_dims
)
  : species_(species)
  , boundary_attribute_(boundary_attribute)
  , num_particles_per_timestep_(num_particles_per_timestep)
  , flux_density_(flux_density)
  , temperature_(temperature)
  , velocity_dims_(velocity_dims)
  , boundary_distribution_(mesh, boundary_attribute)
{}

ParticleContainer MaxwellianParticleBoundarySource::loadParticles(
  RandomNumberGenerator& generator,
  const ParticleOperations& particle_operations,
  double dt
) {
  const double particle_weight =
    flux_density_ * boundary_distribution_.getTotalArea() * dt / num_particles_per_timestep_;

  ParticleContainer particles_on_boundary;
  for (int i = 0; i < num_particles_per_timestep_; i++) {
    const BoundaryPoint boundary_point = boundary_distribution_.generateRandomBoundaryPoint(generator);
    const mfem::Vector velocity = generateVelocity_(boundary_point.outward_unit_normal, generator);
    const mfem::Vector particle_velocity(velocity.GetData(), velocity_dims_);

    particles_on_boundary.addParticle(Particle{
      .position = boundary_point.position,
      .velocity = velocity,
      .element = boundary_point.element,
      .species = species_,
      .weight = particle_weight,
      .is_alive = true,
      .particle_distribution_function_value = distributionFunctionValue(particle_velocity),
    });
  }

  return moveForRandomFractionOfTimestep_(particles_on_boundary, generator, particle_operations, dt);
}

double MaxwellianParticleBoundarySource::distributionFunctionValue(const mfem::Vector& velocity) const {
  const double thermal_speed = std::sqrt(constants::boltzmann_constant * temperature_ / species_.mass);
  const double full_maxwellian_number_density = flux_density_ * std::sqrt(2.0 * M_PI) / thermal_speed;
  const mfem::Vector zero_bulk_velocity({0.0, 0.0, 0.0});
  const mfem::Vector primitive_state = euler::constructPrimitiveState(
    full_maxwellian_number_density, zero_bulk_velocity, temperature_);
  return euler::evaluateMaxwellian(primitive_state, velocity, species_);
}

const Species& MaxwellianParticleBoundarySource::getSpecies() const {
  return species_;
}

int MaxwellianParticleBoundarySource::getBoundaryAttribute() const {
  return boundary_attribute_;
}

int MaxwellianParticleBoundarySource::getNumParticlesPerTimestep() const {
  return num_particles_per_timestep_;
}

double MaxwellianParticleBoundarySource::getFluxDensity() const {
  return flux_density_;
}

double MaxwellianParticleBoundarySource::getTemperature() const {
  return temperature_;
}

mfem::Vector MaxwellianParticleBoundarySource::generateVelocity_(
  const mfem::Vector& outward_unit_normal,
  RandomNumberGenerator& generator
) const {
  mfem::Vector inward_unit_normal = outward_unit_normal;
  inward_unit_normal *= -1.0;
  return generateMaxwellianFluxVelocity(inward_unit_normal, temperature_, species_.mass, generator, velocity_dims_);
}

} // namespace mfpic
