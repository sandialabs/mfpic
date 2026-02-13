#pragma once

#include <libmfpic/ElementFaceContainer.hpp>

#include <memory>

namespace mfem {
  class Vector;
}

namespace mfpic {

struct Particle;

/// Abstract class for applying boundary conditions to particles.
class ParticleBoundary {
public:
  /**
   * @brief Apply the boundary condition to a single particle.
   *
   * @param[in] element_face      Element-local face index upon which particle is incident.
   * @param[in] incoming_particle Incoming particle.
   *
   * @returns Particle with boundary condition applied.
   */
  virtual Particle applyBoundary(int element_face, const Particle& incoming_particle) const = 0;

  /// Dtor.
  virtual ~ParticleBoundary() = default;
};

/// Abstract class for creating particle boundaries.
class ParticleBoundaryFactory {
public:
  /// Common list of parameters passed to all particle boundary factories.
  struct Parameters {
    /// Outward-pointing unit normals for each element-local face on each element.
    std::shared_ptr<ElementFaceContainer<mfem::Vector>> element_face_unit_normal;
  };

  /**
   * @brief Set the boundary attribute to which this particle boundary will apply.
   *
   * @param[in] boundary_attribute Boundary attribute to which this particle boundary will apply.
   */
  void setBoundaryAttribute(int boundary_attribute);

  /**
   * @brief Get the boundary attribute to which this particle boundary will apply.
   *
   * @returns Boundary attribute to which this particle boundary will apply.
   */
  int getBoundaryAttribute() const;

  /**
   * @brief Create the boundary.
   *
   * @param[in] factory_params Parameters used to construct the particle boundary.
   *
   * @returns Particle boundary.
   */
  virtual std::shared_ptr<ParticleBoundary> createBoundary(
    const Parameters& factory_params
  ) const = 0;

private:
  /// Boundary attribute to which this particle boundary will apply.
  int boundary_attribute_;
};

} // namespace mfpic
