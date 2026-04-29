#include <libmfpic/Errors.hpp>
#include <libmfpic/ParticleContainer.hpp>

namespace mfpic {

namespace {

int findSpeciesIndex(const std::vector<Species>& species_list, const Species& species_to_find) {
  const auto species_iterator = std::find(species_list.begin(), species_list.end(), species_to_find);
  if (species_iterator == species_list.end()) {
    errorWithDeveloperMessage("Species not found in list!");
  }
  return std::distance(species_list.begin(), species_iterator);
}

} // namespace

ParticleContainer::ParticleListType::iterator ParticleContainer::begin() {
  return particle_list_.begin();
}

ParticleContainer::ParticleListType::const_iterator ParticleContainer::begin() const {
  return particle_list_.begin();
}

ParticleContainer::ParticleListType::iterator ParticleContainer::end() {
  return particle_list_.end();
}

ParticleContainer::ParticleListType::const_iterator ParticleContainer::end() const {
  return particle_list_.end();
}

int ParticleContainer::numParticles() const {
  return std::ssize(particle_list_);
}

int ParticleContainer::numSpecies() const {
  return std::ssize(species_represented_by_these_particles_);
}

void ParticleContainer::addParticle(Particle particle) {
  particle_list_.push_back(particle);
  const Species& particle_species = particle.species;
  if (std::find(species_represented_by_these_particles_.begin(), species_represented_by_these_particles_.end(), particle_species) == species_represented_by_these_particles_.end()) {
    species_represented_by_these_particles_.push_back(particle_species);
  }
}

void ParticleContainer::addParticles(const ParticleContainer &particles) {
  for (const Particle& particle : particles) {
    addParticle(particle);
  }
}

int ParticleContainer::getParticleSpeciesIndex(const Particle& particle) const {
  const Species& particle_species = particle.species;
  return findSpeciesIndex(species_represented_by_these_particles_, particle_species);
}

void ParticleContainer::cleanOutDeadParticles() {
  std::erase_if(particle_list_, [](const Particle& particle) {return not particle.is_alive;});
}

void ParticleContainer::sortByElementThenSpecies() {
  cleanOutDeadParticles();

  const Particle& particle_with_max_element =
    std::ranges::max(particle_list_, {}, [](const Particle& particle) {return particle.element;});
  const int num_elements = particle_with_max_element.element + 1;

  const int num_species = numSpecies();
  auto flattenElementSpeciesIndex = [num_species](int element, int species_index) {
    return element * num_species + species_index;
  };

  std::vector<int> num_particles_in_element_and_species(num_elements * num_species, 0);
  for (const Particle& particle : particle_list_) {
    num_particles_in_element_and_species[flattenElementSpeciesIndex(
      particle.element,
      getParticleSpeciesIndex(particle)
    )] += 1;
  }

  element_species_bins_ = std::vector<int>(num_elements * num_species + 1, 0);
  std::inclusive_scan(
    num_particles_in_element_and_species.begin(),
    num_particles_in_element_and_species.end(),
    std::next(element_species_bins_.begin())
  );

  std::vector<Particle> sorted_particles(numParticles());
  std::vector<int> new_particle_indices = element_species_bins_;
  for (const Particle& particle : particle_list_) {
    const int new_particle_index = new_particle_indices[flattenElementSpeciesIndex(
      particle.element,
      getParticleSpeciesIndex(particle)
    )]++;
    sorted_particles[new_particle_index] = particle;
  }

  particle_list_ = sorted_particles;
}

std::span<Particle> ParticleContainer::particlesWithElementAndSpecies(int element, const Species& species) {
  sortByElementThenSpecies();

  assert(element < (std::ssize(element_species_bins_) - 1) / numSpecies());

  const int species_index = findSpeciesIndex(species_represented_by_these_particles_, species);
  const int flattened_bin_index = element * numSpecies() + species_index;
  const int element_species_bin_begin = element_species_bins_[flattened_bin_index];
  const int element_species_bin_end = element_species_bins_[flattened_bin_index + 1];
  return std::span(
    std::next(particle_list_.begin(), element_species_bin_begin),
    std::next(particle_list_.begin(), element_species_bin_end)
  );
}

} // namespace mfpic
