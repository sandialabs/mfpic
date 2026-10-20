#pragma once

#include <libmfpic/RandomNumberGenerator.hpp>

namespace mfpic {

class LowFidelityOperations;
class LowFidelityState;
class ParticleContainer;
class ParticleMoments;
class ParticleOperations;

/// Virtual class for performing collisions.
class CollisionOperations {
public:
  /**
   * @brief Apply the collision operator to the particles in place.
   *
   * @param[in]     dt                  Timestep size.
   * @param[in,out] generator           Random number generator.
   * @param[in,out] particles           Particles to which to apply collision operator.
   * @param[in,out] particle_operations Particle operations.
   */
  virtual void performCollisions(
    double dt,
    RandomNumberGenerator& generator,
    ParticleContainer& particles,
    ParticleOperations& particle_operations
  ) const = 0;

  virtual void performCollisions(
    double dt,
    RandomNumberGenerator& generator,
    ParticleContainer& particles,
    ParticleOperations& particle_operations,
    const LowFidelityState& low_fidelity_state,
    const LowFidelityOperations& low_fidelity_operations
  ) const = 0;

  virtual ~CollisionOperations() = default;
};

} // namespace mfpic
