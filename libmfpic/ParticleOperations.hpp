#pragma once

#include <libmfpic/BuildVarianceReductionParametersFromYaml.hpp>
#include <libmfpic/DGEulerOperations.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/ElectromagneticFieldsEvaluator.hpp>
#include <libmfpic/ElementFaceContainer.hpp>
#include <libmfpic/IntegratedCharge.hpp>
#include <libmfpic/LowFidelityOperations.hpp>
#include <libmfpic/LowFidelityState.hpp>
#include <libmfpic/ParticleBoundary.hpp>
#include <libmfpic/ParticleContainer.hpp>

#include <mfem/mfem.hpp>

#include <unordered_map>

#include <unordered_map>

namespace mfpic {

struct ParticleMoments
{
  std::unordered_map<Species, mfem::Vector> number_density;
  std::unordered_map<Species, mfem::DenseMatrix> bulk_velocity;
  std::unordered_map<Species, mfem::Vector> temperature;
};

class ParticleOperations {
public:

  /**
  * @brief Construct a new ParticleOperations
  *
  * @param discretization - Discretization object containing the finite element space
  * @param particle_boundary_factories - List of factories for particle boundaries and attributes to which they apply.
  * @param default_particle_boundary_factory - Factory for particle boundary to apply to uncovered attributes.
  * @param species_map - Map species names to Species structs
  */
  ParticleOperations(
    Discretization &discretization,
    std::vector<std::shared_ptr<ParticleBoundaryFactory>> particle_boundary_factories,
    std::shared_ptr<ParticleBoundaryFactory> default_particle_boundary_factory,
    std::unordered_map<std::string, Species> species_map,
    const int velocity_dims = 3
  );

  ParticleContainer accelerate(
    double dt,
    const ParticleContainer& current_particles,
    const ElectromagneticFieldsEvaluator& field_provider
  ) const;

  /**
   * @brief Move each particle over time @a dt according to its velocity.
   *
   * @param[in] dt                Time interval over which to move particles.
   * @param[in] current_particles Particles before moving.
   *
   * @returns Moved particles.
   */
  ParticleContainer move(
    double dt,
    const ParticleContainer& current_particles
  ) const;

  /**
  * @brief Assemble charges from the particles into the charge density
  *
  * @param current_particles - ParticleContainer of particles
  * @return IntegratedCharge - integrated charge state
  */
  IntegratedCharge assembleCharge(
    const ParticleContainer& current_particles
  ) const;

  /**
  * @brief Assemble charges from the particles into the charge density using a variance reduced formulation
  *
  * The low-fidelity model is the control variate unless Specified LF Moments are set, in which case
  * it is ignored and a Maxwellian at the reference moments is used instead.
  *
  * @param current_particles - ParticleContainer of particles
  * @param low_fidelity_state - \ref LowFidelityState
  * @param low_fidelity_operations - \ref LowFidelityOperations
  * @return IntegratedCharge - integrated charge state
  */
  IntegratedCharge assembleVarianceReducedCharge(
    const ParticleContainer& current_particles,
    const LowFidelityState& low_fidelity_state,
    const LowFidelityOperations& low_fidelity_operations
  ) const;

  /**
  * @brief Assemble variance reduced charge without a low-fidelity model; requires Specified LF Moments
  *
  * @param current_particles - ParticleContainer of particles
  * @return IntegratedCharge - integrated charge state
  */
  IntegratedCharge assembleVarianceReducedCharge(
    const ParticleContainer& current_particles
  ) const;

  /**
   * @brief Compute the particle moments in each element
   *
   * @param[in] particles   \ref ParticleContainer
   *
   * @return Struct of particle number density, bulk velocity, and temperature (species, (element))
   */
  ParticleMoments& getParticleMoments(
    const ParticleContainer& particles
  );

  /**
   * @brief Compute the variance particle moments in each element
   *
   * @param[in] particles   \ref ParticleContainer
   * @param[in] low_fidelity_state \ref LowFidelityState
   * @param[in] low_fidelity_operations \ref LowFidelityOperations
   * @param[in] compute_standard_pic_moments Flag to (re)compute the standard PIC moments
   *
   * @return Struct of variance reduced particle number density, bulk velocity, and temperature (species, (element))
   */
  ParticleMoments& getVarianceReducedParticleMoments(
    const ParticleContainer& particles,
    const LowFidelityState& low_fidelity_state,
    const DGEulerOperations& low_fidelity_operations,
    const bool compute_standard_pic_moments=false
  );

  /**
   * @brief Compute the variance reduced particle moments without a low-fidelity model; requires Specified LF Moments
   *
   * @param[in] particles   \ref ParticleContainer
   * @param[in] compute_standard_pic_moments Flag to (re)compute the standard PIC moments
   *
   * @return Struct of variance reduced particle number density, bulk velocity, and temperature (species, (element))
   */
  ParticleMoments& getVarianceReducedParticleMoments(
    const ParticleContainer& particles,
    const bool compute_standard_pic_moments=false
  );

  /**
   * @brief Get mfem mesh associated with the discretization
   *
   * @return mfem::Mesh
   */
  mfem::Mesh& getMesh() const {return *discretization_.getFeSpace().GetMesh();};

  /**
   * @brief Setter for variance reduction parameters
   */
  void setVarianceReductionParameters(const VarianceReductionParameters& parameters) {
    variance_reduction_parameters_ = parameters;
  }

  /**
   * @brief Getter for variance reduction parameters
   */
  const VarianceReductionParameters getVarianceReductionParameters() const {
    return variance_reduction_parameters_;
  }

  /**
   * @brief Getter for number of velocity dimensions
   */
  int getVelocityDims() const {
    return velocity_dims_;
  }

private:

  /// Control variate integrals over each element, used by the variance reduced moments
  struct ControlVariateMomentIntegrals {
    /// Integral of n, (species, (element))
    std::unordered_map<Species, mfem::Vector> number_density;
    /// Integral of n u, (species, (dimension, element))
    std::unordered_map<Species, mfem::DenseMatrix> momentum;
    /// Integral of n (T + m |u|^2 / (d k_B)), (species, (element))
    std::unordered_map<Species, mfem::Vector> temperature;
  };

  IntegratedCharge assembleVarianceReducedCharge_(
    const ParticleContainer& current_particles,
    const LowFidelityState* low_fidelity_state,
    const LowFidelityOperations* low_fidelity_operations
  ) const;

  ParticleMoments& getVarianceReducedParticleMoments_(
    const ParticleContainer& particles,
    const LowFidelityState* low_fidelity_state,
    const DGEulerOperations* low_fidelity_operations,
    const bool compute_standard_pic_moments
  );

  /**
   * @brief Whether a species has a control variate: listed in Specified LF Moments, or present in the low-fidelity state
   */
  bool hasControlVariate_(
    const Species& species,
    const LowFidelityState* low_fidelity_state
  ) const;

  /**
   * @brief Evaluate the control variate phase-space density at a particle
   */
  double evaluateControlVariate_(
    const Particle& particle,
    const LowFidelityState* low_fidelity_state,
    const LowFidelityOperations* low_fidelity_operations
  ) const;

  /**
   * @brief Flag the elements in which each species uses variance reduction
   *
   * An element is flagged when the species has a control variate and, with Apply Cell Limiting,
   * every particle of that species in the element has |1 - f_LF / f| < 1.
   *
   * @return Map of particle species to flags, (species, (element))
   */
  std::unordered_map<Species, std::vector<char>> computeVarianceReductionElements_(
    const ParticleContainer& particles,
    const LowFidelityState* low_fidelity_state,
    const LowFidelityOperations* low_fidelity_operations
  ) const;

  ControlVariateMomentIntegrals computeControlVariateMomentIntegrals_(
    const LowFidelityState* low_fidelity_state,
    const DGEulerOperations* low_fidelity_operations
  ) const;

  /**
   * @brief Integrated charge of the control variate for one species, per global DOF
   */
  mfem::Vector assembleControlVariateCharge_(
    const Species& species,
    const LowFidelityState* low_fidelity_state,
    const LowFidelityOperations* low_fidelity_operations
  ) const;

  void sumParticleWeights_(
    const ParticleContainer& particles
  ) ;

  /**
   * @brief Compute the number density in each element
   *
   * @param[in] particles   \ref ParticleContainer
   *
   * @return Map of particle species to number density, (species, (element))
   */
  std::unordered_map<Species, mfem::Vector>& getNumberDensity(const ParticleContainer& particles);

  /**
   * @brief Compute the variance reduced number density in each element from the low fidelity state
   *
   * @param[in] particles   \ref ParticleContainer
   * @param[in] low_fidelity_state \ref LowFidelityState
   * @param[in] low_fidelity_operations \ref LowFidelityOperations
   *
   * @return mfem::DenseMatrix of number density for each particle species, (element, species)
   */
    std::unordered_map<Species, mfem::Vector>& getVarianceReducedNumberDensity(
    const ParticleContainer& particles,
    const LowFidelityState* low_fidelity_state,
    const LowFidelityOperations* low_fidelity_operations
  );

  /**
   * @brief Compute the bulk velocity in each element
   *
   * @param[in] particles   \ref ParticleContainer
   * @param[in] sum_weights Optional flag that resums the weights. Default is true.
   *
   * @return Map of particle species to bulk velocity (species, (dimension, element))
   */
  std::unordered_map<Species, mfem::DenseMatrix>& getBulkVelocity(const ParticleContainer& particles, const bool sum_weights = true);

  /**
   * @brief Compute the variance reduced bulk velocity in each element
   *
   * @param[in] particles   \ref ParticleContainer
   * @param[in] low_fidelity_state \ref LowFidelityState
   * @param[in] low_fidelity_operations \ref LowFidelityOperations
   *
   * @return mfem::DenseTensor of bulk velocity for each particle species, (dimension, element, species)
   */
  std::unordered_map<Species, mfem::DenseMatrix>& getVarianceReducedBulkVelocity(
    const ParticleContainer& particles,
    const LowFidelityState* low_fidelity_state,
    const LowFidelityOperations* low_fidelity_operations
  );

  /**
   * @brief Compute the temperature in each element
   *
   * @param[in] particles             \ref ParticleContainer
   * @param[in] sum_weights           Optional flag that resums the weights. Default is true.
   * @param[in] compute_bulk_velocity Optional flag that recomputes the bulk velocity. Default is true.
   *
   * @return Map of particle species to temperature, (species, (element))
   */
  std::unordered_map<Species, mfem::Vector>& getTemperature(const ParticleContainer& particles, const bool sum_weights = true, const bool compute_bulk_velocity = true);

  /**
   * @brief Compute the variance reduced temperature in each element
   *
   * @param[in] particles   \ref ParticleContainer
   * @param[in] low_fidelity_state \ref LowFidelityState
   * @param[in] low_fidelity_operations \ref LowFidelityOperations
   *
   * @return Map of particle species to temperature, (species, (element))
   */
  std::unordered_map<Species, mfem::Vector>& getVarianceReducedTemperature(
    const ParticleContainer& particles,
    const LowFidelityState* low_fidelity_state,
    const LowFidelityOperations* low_fidelity_operations
  );


  /// Discretization object containing the finite element space
  Discretization & discretization_;

  /// Number of faces on each element.
  std::vector<int> num_faces_on_element_;

  /// For each element, for each element-local face, outward-pointing unit normals.
  std::shared_ptr<ElementFaceContainer<mfem::Vector>> element_face_unit_normal_;

  /// For each element, for each element-local face, face centroid dotted with the unit normal.
  ElementFaceContainer<double> element_face_centroid_dot_unit_normal_;

  /// For each element, for each element-local face, the element on the other side of the face.
  ElementFaceContainer<int> element_face_other_element_;

  /// Particle boundaries.
  ElementFaceContainer<std::shared_ptr<ParticleBoundary>> particle_boundaries_;

  /// Particle moments
  ParticleMoments particle_moments_;

  /// Particle moments
  ParticleMoments variance_reduced_particle_moments_;

  /// Variance reduction parameters
  VarianceReductionParameters variance_reduction_parameters_;

  /// Particle sum of weights
  std::unordered_map<Species, mfem::Vector> sum_of_weights_;

  /// Mesh dimension
  const int dim_;

  /// Velocity dimension
  int velocity_dims_ = 3;

  /// Elements in which each species uses variance reduction, set by getVarianceReducedParticleMoments
  std::unordered_map<Species, std::vector<char>> variance_reduction_elements_;

  /// Control variate moment integrals, set by getVarianceReducedParticleMoments
  ControlVariateMomentIntegrals control_variate_moment_integrals_;
};

} // namespace mfpic
