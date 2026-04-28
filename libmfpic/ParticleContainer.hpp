#pragma once

#include <libmfpic/Particle.hpp>

namespace mfpic {

class ParticleContainer {
public:
  using ParticleListType = std::vector<Particle>;

  ParticleListType::iterator begin();

  ParticleListType::const_iterator begin() const;

  ParticleListType::iterator end();

  ParticleListType::const_iterator end() const;

  int numParticles() const;

  void addParticle(Particle particle);

  void addParticles(const ParticleContainer& particles);

  int getParticleSpeciesIndex(const Particle& particle) const;

  void cleanOutDeadParticles();

private:
  std::vector<Particle> particle_list_;

  std::vector<Species> species_represented_by_these_particles_;
};

} // namespace mfpic
