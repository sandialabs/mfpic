#pragma once

#include <libmfpic/ParticleContainer.hpp>
#include <libmfpic/ParticleOperations.hpp>
#include <libmfpic/RandomNumberGenerator.hpp>

namespace mfpic {

/// Abstract class for injecting particles through a mesh boundary each time step.
class ParticleBoundarySource {
public:
  /**
   * @brief Create the particles injected through the boundary over one time step.
   *
   * @param[in,out] generator           Random number generator.
   * @param[in]     particle_operations Particle operations used to move the injected particles.
   * @param[in]     dt                  Time step size.
   *
   * @returns Container of injected particles.
   */
  virtual ParticleContainer loadParticles(
    RandomNumberGenerator& generator,
    const ParticleOperations& particle_operations,
    double dt
  ) = 0;

  /// Dtor.
  virtual ~ParticleBoundarySource();

protected:
  /**
   * @brief Move particles on the boundary inward as if each crossed it at a uniformly random time during the step.
   *
   * Each particle is moved for a uniformly random fraction of @a dt.
   * The electric field is not applied during this partial step.
   * Particles that leave the domain during the partial step through an absorbing boundary are removed.
   *
   * @param[in]     particles_on_boundary Particles located on the boundary.
   * @param[in,out] generator             Random number generator.
   * @param[in]     particle_operations   Particle operations used to move the particles.
   * @param[in]     dt                    Time step size.
   *
   * @returns Moved particles that are still alive.
   */
  static ParticleContainer moveForRandomFractionOfTimestep_(
    const ParticleContainer& particles_on_boundary,
    RandomNumberGenerator& generator,
    const ParticleOperations& particle_operations,
    double dt
  );
};

} // namespace mfpic
