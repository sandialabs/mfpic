#pragma once

#include <libmfpic/BoundaryDistribution.hpp>
#include <libmfpic/ParticleBoundarySource.hpp>
#include <libmfpic/Species.hpp>

#include <mfem/mfem.hpp>

#include <memory>

namespace mfpic {

/**
 * @brief Boundary source injecting particles that cross the boundary from a non-drifting Maxwellian gas.
 *
 * Particles are spread uniformly over the boundary faces.
 * Normal velocities follow the flux-weighted distribution \f$ f(v_n) \propto v_n \exp(-v_n^2 / 2\sigma^2) \f$
 * and tangential velocities are Maxwellian, so that the injected particles form a half-Maxwellian.
 * Next to the boundary its number density is \f$ \Gamma \sqrt{\pi/2} / \sigma \f$ with \f$ \sigma = \sqrt{k_B T / m} \f$.
 */
class MaxwellianParticleBoundarySource : public ParticleBoundarySource {
public:
  /**
   * @brief Ctor.
   *
   * @throws If @a mesh contains no boundary faces with @a boundary_attribute.
   *
   * @param[in] mesh                       Mesh through whose boundary particles are injected.
   * @param[in] boundary_attribute         Boundary attribute of the faces through which particles are injected.
   * @param[in] species                    Species of the injected particles.
   * @param[in] num_particles_per_timestep Number of macroparticles injected per time step.
   * @param[in] flux_density               Number of physical particles crossing the boundary per unit area per unit time.
   * @param[in] temperature                Temperature of the Maxwellian gas.
   * @param[in] velocity_dims              Number of velocity dimensions.
   */
  MaxwellianParticleBoundarySource(
    std::shared_ptr<mfem::Mesh> mesh,
    int boundary_attribute,
    const Species& species,
    int num_particles_per_timestep,
    double flux_density,
    double temperature,
    int velocity_dims
  );

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
  ) override;

  /**
   * @brief Evaluate the phase-space density of the injected particles, used for variance reduction.
   *
   * @param[in] velocity Particle velocity, of length equal to the number of velocity dimensions.
   *
   * @returns Value of the phase-space density at @a velocity next to the boundary.
   */
  double distributionFunctionValue(const mfem::Vector& velocity) const;

  /// Get the species of the injected particles.
  const Species& getSpecies() const;

  /// Get the boundary attribute of the faces through which particles are injected.
  int getBoundaryAttribute() const;

  /// Get the number of macroparticles injected per time step.
  int getNumParticlesPerTimestep() const;

  /// Get the number of physical particles crossing the boundary per unit area per unit time.
  double getFluxDensity() const;

  /// Get the temperature of the Maxwellian gas.
  double getTemperature() const;

private:
  /**
   * @brief Draw the velocity of a particle entering the domain through the boundary.
   *
   * @param[in]     outward_unit_normal Unit normal of the boundary pointing out of the domain, of length 3.
   * @param[in,out] generator           Random number generator.
   *
   * @returns Velocity of length 3 with a negative component along @a outward_unit_normal.
   */
  mfem::Vector generateVelocity_(const mfem::Vector& outward_unit_normal, RandomNumberGenerator& generator) const;

  /// Species of the injected particles.
  Species species_;

  /// Boundary attribute of the faces through which particles are injected.
  int boundary_attribute_;

  /// Number of macroparticles injected per time step.
  int num_particles_per_timestep_;

  /// Number of physical particles crossing the boundary per unit area per unit time.
  double flux_density_;

  /// Temperature of the Maxwellian gas.
  double temperature_;

  /// Number of velocity dimensions.
  int velocity_dims_;

  /// Distribution of points on the faces through which particles are injected.
  BoundaryDistribution boundary_distribution_;
};

} // namespace mfpic
