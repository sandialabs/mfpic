#pragma once

#include <libmfpic/Discretization.hpp>
#include <libmfpic/ElectromagneticFieldsEvaluator.hpp>
#include <libmfpic/ElementFaceContainer.hpp>
#include <libmfpic/IntegratedCharge.hpp>
#include <libmfpic/ParticleBoundary.hpp>
#include <libmfpic/ParticleContainer.hpp>

#include <mfem/mfem.hpp>

namespace mfpic {

class ParticleOperations {
public:

  /**
  * @brief Construct a new ParticleOperations 
  *
  * @param discretization - Discretization object containing the finite element space 
  * @param particle_boundary_factories - List of factories for particle boundaries and attributes to which they apply.
  * @param default_particle_boundary_factory - Factory for particle boundary to apply to uncovered attributes.
  */
  ParticleOperations(
    Discretization &discretization,
    std::vector<std::shared_ptr<ParticleBoundaryFactory>> particle_boundary_factories,
    std::shared_ptr<ParticleBoundaryFactory> default_particle_boundary_factory
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

private: 

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

  /// Mesh dimension
  const int dim_;

};

} // namespace mfpic
