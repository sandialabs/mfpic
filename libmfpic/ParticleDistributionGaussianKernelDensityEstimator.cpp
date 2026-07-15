#include <libmfpic/Constants.hpp>
#include <libmfpic/ParticleDistributionGaussianKernelDensityEstimator.hpp>

#include <ranges>

namespace mfpic {

// LTM do this better
void ParticleDistributionGaussianKernelDensityEstimator::setUpEstimation(
  ParticleContainer& particles,
  ParticleOperations& particle_operations
) {
  particles.sortByElementThenSpecies();

  number_densities_ = particle_operations.getNumberDensity(particles);

  const std::unordered_map<Species, mfem::Vector>& temperatures =
    particle_operations.getTemperature(particles, false, true);

  bandwidths_ = number_densities_;
  std::unordered_map<Species, mfem::Vector> sum_of_weights = number_densities_;
  for (mfem::Vector& weight_sum : std::views::values(sum_of_weights)) {
    weight_sum = 0.0;
  }
  std::unordered_map<Species, mfem::Vector> sum_of_squared_weights = number_densities_;
  for (mfem::Vector& weight_sum : std::views::values(sum_of_squared_weights)) {
    weight_sum = 0.0;
  }
  for (const Particle& particle : particles) {
    const double particle_weight = particle.weight;
    sum_of_weights.at(particle.species)(particle.element) += particle_weight;
    sum_of_squared_weights.at(particle.species)(particle.element) += particle_weight * particle_weight;
  }

  for (auto& [species, bandwidths] : bandwidths_) {
    for (int element = 0; element < bandwidths.Size(); element++) {
      const double effective_particle_count =
        std::pow(sum_of_weights.at(species)(element), 2.0) / sum_of_squared_weights.at(species)(element);
      constexpr int dimensions = 3;
      const double sigma = std::sqrt(constants::boltzmann_constant * temperatures.at(species)(element) / species.mass);
      bandwidths_.at(species)(element) = sigma * std::pow(effective_particle_count, -1.0 / (dimensions + 4));
    }
  }
}

double ParticleDistributionGaussianKernelDensityEstimator::estimateDensity(
  const mfem::Vector&,
  int element,
  const mfem::Vector& velocity,
  const Species& species,
  const ParticleContainer& particles,
  const ParticleOperations& particle_operations
) const {
  assert(particles.isSortedByElementThenSpecies());

  const double number_density_in_element = number_densities_.at(species)(element);
  mfem::Mesh& mesh = particle_operations.getMesh();

  const double element_volume = mesh.GetElementVolume(element);
  double density_estimate = 0.0;
  for (const Particle& particle : particles.particlesWithElementAndSpecies(element, species)) {
    const double bandwidth = bandwidths_.at(species)(element);
    mfem::Vector relative_velocity = velocity;
    relative_velocity -= particle.velocity;
    density_estimate +=
      particle.weight / element_volume / number_density_in_element /
      std::pow(2.0 * M_PI * bandwidth * bandwidth, 3./2.) *
      std::exp(- (relative_velocity*relative_velocity) / 2.0 / bandwidth / bandwidth);
  }

  return density_estimate;
}

} // namespace mfpic
