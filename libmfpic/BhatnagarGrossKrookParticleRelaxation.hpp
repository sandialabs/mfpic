#pragma once

#include <libmfpic/CollisionOperations.hpp>
#include <libmfpic/Species.hpp>

namespace mfpic {

class BhatnagarGrossKrookParticleRelaxation : public CollisionOperations {
public:
  BhatnagarGrossKrookParticleRelaxation(const Species& species, double relaxation_frequency);

  virtual void performCollisions(
    RandomNumberGenerator& generator,
    ParticleContainer& particles,
    LowFidelityState&
  ) const override;

private:
  const Species species_;

  const double relaxation_frequency_;
};

} // namespace mfpic
