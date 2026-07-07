#pragma once

namespace mfpic {

class ParticleContainer;
class ParticleOperations;
class RandomNumberGenerator;

class CollisionOperations {
public:
  virtual void performCollisions(
    double dt,
    RandomNumberGenerator& generator,
    ParticleContainer& particles,
    ParticleOperations& particle_operations
  ) const = 0;
};

} // namespace mfpic
