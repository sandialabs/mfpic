#pragma once

#include <libmfpic/CollisionOperations.hpp>

namespace mfpic {

class BGKRelaxationOperations : public CollisionOperations {
public:
  BGKRelaxationOperations(double collision_frequency);

  virtual void performCollisions(
    RandomNumberGenerator& generator,
    ParticleContainer& particles,
    LowFidelityState&
  ) const override;

private:
  const double collision_frequency_;
};

} // namespace mfpic
