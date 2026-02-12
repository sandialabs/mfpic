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

private:
  std::vector<Particle> particle_list_;
};

} // namespace mfpic
