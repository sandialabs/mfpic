#pragma once

#include <libmfpic/CollisionOperations.hpp>
#include <libmfpic/Particle.hpp>
#include <libmfpic/Species.hpp>

#include <vector>

namespace mfpic {

/**
 * @brief Operations for the Bhatnagar-Gross-Krook collision operator.
 *
 * This operator takes the form \f$ \frac{df}{dt} = \nu (f - f_0) \f$,
 * where \f$ f_0 \f$ is the equilibrium distribution
 * and \f$ \nu \f$ is the collision frequency.
 * The equilibrium distribution in each cell is a drifting Maxwellian
 * with the same number density, bulk velocity, and temperature
 * as the macroparticles in the cell.
 *
 * In each timestep, each macroparticle relaxes independently with probability
 * \f$ 1 - e^{-\nu \Delta t} \f$. Relaxing particles have their velocities
 * resampled from the cell Maxwellian. Mass, momentum, and energy are conserved
 * in expectation.
 */
class BGKRelaxationOperations : public CollisionOperations {
public:
  /**
   * @brief Ctor.
   *
   * @param[in] collision_frequency Collision frequency \f$ \nu \f$.
   * @param[in] species_to_relax    Species to which to apply the BGK operator.
   * @param[in] averaging_timesteps e-folding time, in timesteps, of the exponential moving average of the
   *                                equilibrium moments (1 = no time averaging).
   * @param[in] smoothing_length    Elements whose centers are within this distance are averaged together
   *                                (0 = no spatial averaging).
   * @param[in] smoothing_order     Polynomial order of the local least-squares fit over the smoothing window
   *                                (0 = flat average, 1 = linear, 2 = quadratic).
   * @param[in] use_vr_moments      Draw relaxed velocities and stamp f_HF from the variance-reduced particle
   *                                moments (energy-conserving T and u) when a matching fluid species exists.
   *                                Smoothing options are ignored in this case.
   * @param[in] stamp_raw_density   With use_vr_moments, stamp f_HF with the element's raw particle number density
   *                                (which BGK conserves exactly) instead of the variance-reduced one.
   * @param[in] ratio_stamp         Update f_HF of a relaxed particle multiplicatively, f_HF <- f_HF M(v_new) / M(v_old),
   *                                with M the stamping Maxwellian, so the particle keeps its phase-space volume
   *                                instead of being renormalized to the element's current density.
   * @param[in] density_scale_stride If positive, stamp relaxed particles with f_HF = c M(v_new) / n, where the
   *                                per-element density scale c is the mean equilibrium density over the previous
   *                                block of this many calls (held fixed within a block) instead of the current one.
   * @param[in] invariant_constraint After each call, rescale f_HF of every particle in the element by a maximum
   *                                cross-entropy factor exp(-lambda . R(v)), R = {1, xi, |xi|^2}, so that the
   *                                equilibrium masses w M(v) / f_HF carry the same R moments as before the call
   *                                (Sadr & Hadjiconstantinou 2023). Relaxed particles start from the mean
   *                                pre-call equilibrium mass per weight of the relaxed particles, or from their
   *                                KDE mass when kde_neighbors > 0.
   * @param[in] kde_neighbors       If positive, stamp a relaxed particle with f_HF = M(v_new) / W, where W is the
   *                                weight-averaged local-basis weight M(v_j) / f_HF_j of the kde_neighbors
   *                                pre-call particles in the element nearest to v_new in thermal-speed units
   *                                (Al-Mohssen & Hadjiconstantinou 2010). Takes precedence over ratio_stamp.
   * @param[in] relaxed_mass_renormalization Rescale the f_HF of an element's relaxed particles by a common factor so
   *                                their total local-basis equilibrium mass w M(v) / f_HF equals its pre-call value.
   * @param[in] kde_unit_scale      Rescale the f_HF of an element's relaxed particles by a common factor so their
   *                                weight-averaged local-basis weight M(v) / f_HF is 1, i.e. take the shape from the
   *                                KDE and the scale from the stamping Maxwellian. Overrides relaxed_mass_renormalization.
   * @param[in] log_cell_energy     Print per-element BGK energy change, kinetic energy, VR-active flag, and
   *                                maximum |1 - f_LF/f_HF| on every call.
   */
  BGKRelaxationOperations(double collision_frequency, const Species& species_to_relax,
    double averaging_timesteps = 1.0, double smoothing_length = 0.0, int smoothing_order = 0,
    bool use_vr_moments = false, bool log_cell_energy = false, bool stamp_raw_density = false,
    bool ratio_stamp = false, int density_scale_stride = 0, bool invariant_constraint = false,
    int kde_neighbors = 0, bool relaxed_mass_renormalization = false,
    bool kde_unit_scale = false);

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

  virtual void performCollisions(
    double dt,
    RandomNumberGenerator& generator,
    ParticleContainer& particles,
    ParticleOperations& particle_operations,
    const LowFidelityState& low_fidelity_state,
    const LowFidelityOperations& low_fidelity_operations
  ) const override;

private:
  std::vector<std::vector<Particle*>> selectRelaxingParticlesPerElement_(
    double dt,
    RandomNumberGenerator& generator,
    ParticleContainer& particles,
    int num_elements
  ) const;

  /// Collision frequency \f$ \nu \f$.
  const double collision_frequency_;

  /// Species to which to apply the BGK operator.
  const Species species_to_relax_;

  const double averaging_timesteps_;

  const double smoothing_length_;

  const int smoothing_order_;

  const bool use_vr_moments_;

  const bool log_cell_energy_;

  const bool stamp_raw_density_;

  const bool ratio_stamp_;

  const int density_scale_stride_;

  const bool invariant_constraint_;

  const int kde_neighbors_;

  const bool relaxed_mass_renormalization_;

  const bool kde_unit_scale_;

  mutable long long relaxed_mass_renormalization_failures_ = 0;

  /// Weight-averaged local-basis weight of the kde_neighbors_ samples nearest to normalized velocity xi.
  double kdeLocalWeight_(const mfem::Vector& xi, const std::vector<mfem::Vector>& velocities,
    const std::vector<double>& local_weights, const std::vector<double>& particle_weights) const;

  /// Particles of the relaxed species in each element at the last call (filled only with invariant_constraint_ or KDE).
  mutable std::vector<std::vector<Particle*>> particles_per_element_;

  mutable long long invariant_constraint_failures_ = 0;
  mutable int invariant_constraint_max_iterations_ = 0;

  /// Solve the per-element cross-entropy problem for the five collisional invariants and restamp f_HF.
  void constrainInvariants_(int element, const mfem::Vector& stamp_state, int velocity_dims,
    const std::vector<double>& equilibrium_masses_before, const std::vector<mfem::Vector>& velocities_before,
    const std::vector<Particle*>& relaxed_particles) const;

  /// Per-element density used to stamp f_HF when density_scale_stride_ > 0, and its running block sum.
  mutable mfem::Vector density_scale_;
  mutable mfem::Vector density_scale_sum_;
  mutable int density_scale_calls_ = 0;

  /// Add this call's equilibrium densities to the block sum and refresh density_scale_ at the end of each block.
  void updateDensityScale_(const mfem::Vector& number_densities) const;

  /// Number of ratio-stamp updates that fell back to n M(v_new) because the ratio was not finite and positive.
  mutable long long ratio_stamp_fallbacks_ = 0;

  /// Per-element kinetic energy change of the last relaxToEquilibria_ call.
  mutable std::vector<double> cell_energy_change_;

  /// Relax particles in each element toward the Maxwellian with that element's primitive state
  /// (n, u_x, u_y, u_z, T), and stamp their f_HF from it. Elements with an empty state are skipped.
  void relaxToEquilibria_(
    double dt,
    RandomNumberGenerator& generator,
    ParticleContainer& particles,
    int velocity_dims,
    const std::vector<mfem::Vector>& primitive_states,
    const mfem::Vector* stamp_number_densities = nullptr
  ) const;

  /// Per-element conserved densities: number, momentum (3), and n<|v|^2>.
  struct ConservedDensities {
    mfem::Vector number;
    mfem::DenseMatrix momentum;
    mfem::Vector speed_squared;
  };

  /// Space-smoothed conserved densities from the current particles.
  ConservedDensities smoothedConservedDensities_(const ParticleMoments& moments, const mfem::Mesh& mesh, int velocity_dims) const;

  /// Time-averaged, space-smoothed conserved densities defining the equilibrium.
  mutable ConservedDensities averaged_densities_;

  mutable bool averaged_densities_initialized_ = false;

  /// For each element, the elements (including itself) within the smoothing length.
  mutable std::vector<std::vector<int>> smoothing_neighbors_;
};

} // namespace mfpic
