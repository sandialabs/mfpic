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
}

void ParticleContainer::addParticles(const ParticleContainer &particles) {
  for (const Particle& particle : particles) {
    addParticle(particle);
  }
}

} // namespace mfpic
