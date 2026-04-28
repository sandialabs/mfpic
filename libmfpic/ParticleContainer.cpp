#include <libmfpic/ParticleContainer.hpp>

namespace mfpic {

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
  const auto species_iterator = std::find(species_represented_by_these_particles_.begin(), species_represented_by_these_particles_.end(), particle_species);
  return std::distance(species_represented_by_these_particles_.begin(), species_iterator);
}

void ParticleContainer::cleanOutDeadParticles() {
  std::erase_if(particle_list_, [](const Particle& particle) {return not particle.is_alive;});
}

} // namespace mfpic
