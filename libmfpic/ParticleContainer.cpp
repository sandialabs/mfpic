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
  is_sorted_ = false;
  return particle_list_.begin();
}

ParticleContainer::ParticleListType::const_iterator ParticleContainer::begin() const {
  return particle_list_.begin();
}

ParticleContainer::ParticleListType::iterator ParticleContainer::end() {
  is_sorted_ = false;
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
  is_sorted_ = false;
}

void ParticleContainer::addParticles(const ParticleContainer &particles) {
  for (const Particle& particle : particles) {
    addParticle(particle);
  }
  is_sorted_ = false;
}

int ParticleContainer::getParticleSpeciesIndex(const Particle& particle) const {
  const Species& particle_species = particle.species;
  return findSpeciesIndex(species_represented_by_these_particles_, particle_species);
}

void ParticleContainer::cleanOutDeadParticles() {
  std::erase_if(particle_list_, [](const Particle& particle) {return not particle.is_alive;});
}

void ParticleContainer::sortByElementThenSpecies() {
  if (is_sorted_) return;

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

  element_species_bin_offsets_ = std::vector<int>(num_elements * num_species);
  std::exclusive_scan(
    num_particles_in_element_and_species.begin(),
    num_particles_in_element_and_species.end(),
    element_species_bin_offsets_.begin(),
    0
  );

  std::vector<Particle> sorted_particles(numParticles());
  std::vector<int> new_particle_indices = element_species_bin_offsets_;
  for (const Particle& particle : particle_list_) {
    const int new_particle_index = new_particle_indices[flattenElementSpeciesIndex(
      particle.element,
      getParticleSpeciesIndex(particle)
    )]++;
    sorted_particles[new_particle_index] = particle;
  }

  particle_list_ = sorted_particles;
  is_sorted_ = true;
}

} // namespace mfpic
