#pragma once

#include <libmfpic/ElementFaceContainer.hpp>
#include <libmfpic/ParticleBoundary.hpp>

#include <mfem/mfem.hpp>

#include <forward_list>

namespace mfpic {

/// Particle boundary that handles translating particles across a periodic mesh.
class PeriodicParticleBoundary : public ParticleBoundary {
public:
  /**
   * @brief Create a map from periodic faces to this periodic particle boundary handler.
   *
   * Periodic boundaries do not use a factory like other particle boundaries because they are not requested by users.
   * Rather, periodic boundaries are a property of the mesh, and so it should always be checked whether they are needed.
   *
   * @param[in] mesh Mesh to check for periodic boundaries.
   *
   * @returns A map from periodic faces to periodic particle boundaries. If no such faces exist, this will be empty.
   */
  static ElementFaceContainer<std::shared_ptr<ParticleBoundary>> generatePeriodicParticleBoundaries(mfem::Mesh& mesh);

  /**
   * @brief Translate a particle across a periodic mesh and update its element.
   *
   * @param[in] element_face      Element-local face index upon which particle is incident.
   * @param[in] incoming_particle Incoming particle.
   *
   * @returns Translated particle.
   */
  virtual Particle applyBoundary(int element_face, const Particle& incoming_particle) const override;

private:
  /// Collections of transformation objects for each of the periodic faces.
  ElementFaceContainer<mfem::FaceElementTransformations> face_element_transformations_;

  /// A dump of element transformation objects to which the collections in @ref face_element_transformations_ point.
  std::forward_list<mfem::IsoparametricTransformation> element_transformation_dump_;
};

} // namespace mfpic
