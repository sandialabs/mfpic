#include <libmfpic/ParticleContainer.hpp>

namespace mfpic {

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
  const auto species_iterator = std::find(species_represented_by_these_particles_.begin(), species_represented_by_these_particles_.end(), particle_species);
  return std::distance(species_represented_by_these_particles_.begin(), species_iterator);
}

void ParticleContainer::cleanOutDeadParticles() {
  std::erase_if(particle_list_, [](const Particle& particle) {return not particle.is_alive;});
}

} // namespace mfpic
