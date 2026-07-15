#pragma once

#include <libmfpic/ParticleDistributionDensityEstimator.hpp>

namespace mfpic {

class ParticleDistributionGaussianKernelDensityEstimator : public ParticleDistributionDensityEstimator {
public:
  virtual void setUpEstimation(
    ParticleContainer& particles,
    ParticleOperations& particle_operations
  ) override;

  virtual double estimateDensity(
    const mfem::Vector& x,
    int element,
    const mfem::Vector& velocity,
    const Species& species,
    const ParticleContainer& particles,
    const ParticleOperations& particle_operations
  ) const override;

private:
  std::unordered_map<Species, mfem::Vector> bandwidths_;

  std::unordered_map<Species, mfem::Vector> number_densities_;
};

} // namespace mfpic
