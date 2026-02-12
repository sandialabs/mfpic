#pragma once

#include <libmfpic/ParticleBoundary.hpp>

namespace mfpic {

/**
 * @brief Specularly reflects particles incident on faces.
 */
class ReflectingParticleBoundary : public ParticleBoundary {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] element_face_unit_normal Outward-pointing unit normals for each element and face.
   */
  ReflectingParticleBoundary(std::shared_ptr<ElementFaceContainer<mfem::Vector>> element_face_unit_normal);

  /**
   * @brief Reflect a particle.
   *
   * @param[in] element_face      Element-local face index upon which particle is incident.
   * @param[in] incoming_particle Incoming particle.
   *
   * @returns Reflected particle.
   */
  virtual Particle applyBoundary(int element_face, const Particle &incoming_particle) const override;

private:
  /// Outward-pointing unit normals for each element and face.
  std::shared_ptr<ElementFaceContainer<mfem::Vector>> element_face_unit_normal_;
};

/// Factory used to create a reflecting particle boundary.
class ReflectingParticleBoundaryFactory : public ParticleBoundaryFactory {
public:
  /**
   * @brief Create the boundary.
   *
   * @param[in] factory_params Parameters used to construct the particle boundary, unit normals in this case.
   *
   * @returns Particle boundary.
   */
  virtual std::shared_ptr<ParticleBoundary> createBoundary(
    const Parameters& factory_params
  ) const override;
};

} // namespace mfpic
