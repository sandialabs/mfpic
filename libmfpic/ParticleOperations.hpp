#pragma once

#include <libmfpic/DGEulerOperations.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/ElectromagneticFieldsEvaluator.hpp>
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
    std::unordered_map<std::string, Species> species_map
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
  * @param current_particles - ParticleContainer of particles
  * @return IntegratedCharge - integrated charge state
  */
  IntegratedCharge assembleVarianceReducedCharge(
    const ParticleContainer& current_particles,
    const LowFidelityState& low_fidelity_state,
    const LowFidelityOperations& low_fidelity_operations
  ) const;

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
   *
   * @return mfem::DenseMatrix of number density for each particle species, (element, species)
   */
  mfem::DenseMatrix& getVarianceReducedNumberDensity(
    const ParticleContainer& particles,
    const LowFidelityState& low_fidelity_state,
    const DGEulerOperations& low_fidelity_operations
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
   * @brief Compute the bulk velocity in each element
   *
   * @param[in] particles   \ref ParticleContainer
   * @param[in] sum_weights Optional flag that resums the weights. Default is true.
   *
   * @return mfem::DenseTensor of bulk velocity for each particle species, (dimension, element, species)
   */
  mfem::DenseTensor& getVarianceReducedBulkVelocity(
    const ParticleContainer& particles, 
    const LowFidelityState& low_fidelity_state,
    const DGEulerOperations& low_fidelity_operations
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
   * @brief Get mfem mesh associated with the discretization
   *
   * @return mfem::Mesh 
   */
  mfem::Mesh& getMesh() const {return *discretization_.getFeSpace().GetMesh();};

private: 

  void sumParticleWeights_(
    const ParticleContainer& particles
  ) ;

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

  /// Particle number density
  std::unordered_map<Species, mfem::Vector> particle_number_density_;

  /// Particle bulk velocity
  std::unordered_map<Species, mfem::DenseMatrix> particle_bulk_velocity_;

  /// Particle temperature
  std::unordered_map<Species, mfem::Vector> particle_temperature_;

  /// Variance reduced particle number density
  mfem::DenseMatrix variance_reduced_particle_number_density_;

  /// Variance reduced particle bulk velocity
  mfem::DenseTensor variance_reduced_particle_bulk_velocity_;

  /// Variance reduced particle temperature
  mfem::DenseMatrix variance_reduced_particle_temperature_;

  /// Particle sum of weights
  std::unordered_map<Species, mfem::Vector> sum_of_weights_;

  /// Mesh dimension
  const int dim_;

};

} // namespace mfpic
