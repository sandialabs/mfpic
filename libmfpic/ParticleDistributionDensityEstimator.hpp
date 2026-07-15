#pragma once

#include <libmfpic/ParticleOperations.hpp>

namespace mfpic {

class ParticleDistributionDensityEstimator {
public:
  virtual void setUpEstimation(
    ParticleContainer& particles,
    ParticleOperations& particle_operations
  ) = 0;

  virtual double estimateDensity(
    const mfem::Vector& x,
    int element,
    const mfem::Vector& velocity,
    const Species& species,
    const ParticleContainer& particles,
    const ParticleOperations& particle_operations
  ) const = 0;
};

} // namespace mfpic
