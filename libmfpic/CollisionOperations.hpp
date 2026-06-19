#pragma once

namespace mfpic {

class ParticleContainer;
class LowFidelityState;

class CollisionOperations {
public:
  virtual void performCollisions(
    RandomNumberGenerator& generator,
    ParticleContainer& particles,
    LowFidelityState& low_fidelity_state
  ) const = 0;
};

} // namespace mfpic
