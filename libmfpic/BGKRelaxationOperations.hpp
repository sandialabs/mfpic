#pragma once

#include <libmfpic/CollisionOperations.hpp>
#include <libmfpic/Species.hpp>

namespace mfpic {

class BGKRelaxationOperations : public CollisionOperations {
public:
  BGKRelaxationOperations(double collision_frequency, const Species& species_to_relax);

  virtual void performCollisions(
    double dt,
    RandomNumberGenerator& generator,
    ParticleContainer& particles,
    ParticleOperations& particle_operations
  ) const override;

private:
  const double collision_frequency_;

  const Species species_to_relax_;
};

} // namespace mfpic
