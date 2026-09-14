#pragma once

#include <libmfpic/CollisionOperations.hpp>
#include <libmfpic/Species.hpp>

namespace mfpic {

/**
 * @brief Operations for the Bhatnagar-Gross-Krook collision operator.
 *
 * This operator takes the form \f$ \frac{df}{dt} = \nu (f - f_0) \f$,
 * where \f$ f_0 \f$ is the equilibrium distribution
 * and \f$ \nu \f$ is the collision frequency.
 * We assume that the equilibrium distribution in each cell is a drifting Maxwellian
 * with the same number density, bulk velocity, and temperature
 * as the extant macroparticles in the cell.
 */
class BGKRelaxationOperations : public CollisionOperations {
public:
  /**
   * @brief Ctor.
   *
   * @param[in] collision_frequency Collision frequency \f$ \nu \f$.
   * @param[in] species_to_relax    Species to which to apply the BGK operator.
   */
  BGKRelaxationOperations(double collision_frequency, const Species& species_to_relax);

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
  ) const override;

private:
  /// Collision frequency \f$ \nu \f$.
  const double collision_frequency_;

  /// Species to which to apply the BGK operator.
  const Species species_to_relax_;
};

} // namespace mfpic
