#include <libmfpic/BGKRelaxationOperations.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/GenerateMaxwellianVelocity.hpp>
#include <libmfpic/ParticleOperations.hpp>

namespace mfpic {

BGKRelaxationOperations::BGKRelaxationOperations(double collision_frequency, const Species& species_to_relax)
: collision_frequency_(collision_frequency),
  species_to_relax_(species_to_relax)
{
}

void BGKRelaxationOperations::performCollisions(
  double dt,
  RandomNumberGenerator& generator,
  ParticleContainer& particles,
  ParticleOperations& particle_operations
) const {
  const std::unordered_map<Species, mfem::Vector>& number_densities = particle_operations.getNumberDensity(particles);
  std::unordered_map<Species, mfem::DenseMatrix>& bulk_velocities = particle_operations.getBulkVelocity(particles);
  constexpr bool recompute_lower_order_moments = false;
  const std::unordered_map<Species, mfem::Vector>& temperatures = particle_operations.getTemperature(
    particles,
    recompute_lower_order_moments,
    recompute_lower_order_moments
  );

  std::uniform_real_distribution<double> uniform_unit_interval_distribution;

  for (Particle& particle : particles) {
    if (not (particle.is_alive and particle.species == species_to_relax_)) {
      continue;
    }

    const double relaxation_probability = -std::expm1(-collision_frequency_ * dt);
    if (uniform_unit_interval_distribution(generator) < relaxation_probability) {
      const int element = particle.element;
      const mfem::Vector bulk_velocity(bulk_velocities.at(species_to_relax_).GetColumn(element), 3);
      const double temperature = temperatures.at(species_to_relax_)[element];
      particle.velocity = generateMaxwellianVelocity(bulk_velocity, temperature, species_to_relax_.mass, generator);

      const mfem::Vector primitive_state {
        number_densities.at(species_to_relax_)(element),
        bulk_velocity(0),
        bulk_velocity(1),
        bulk_velocity(2),
        temperature};

      particle.particle_distribution_function_value = euler::evaluateMaxwellian(primitive_state, particle.velocity, species_to_relax_);
    }
  }
}

} // namespace mfpic
