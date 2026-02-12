#include <libmfpic/ElementFaceContainer.hpp>
#include <libmfpic/Particle.hpp>
#include <libmfpic/ReflectingParticleBoundary.hpp>

namespace mfpic {

ReflectingParticleBoundary::ReflectingParticleBoundary(
  std::shared_ptr<ElementFaceContainer<mfem::Vector>> element_face_unit_normal
) :
  element_face_unit_normal_(element_face_unit_normal)
{
}

Particle ReflectingParticleBoundary::applyBoundary(int element_face, const Particle &incoming_particle) const {
  Particle reflected_particle = incoming_particle;

  reflected_particle.velocity.Add(
    -2.0 * (reflected_particle.velocity * element_face_unit_normal_->at(reflected_particle.element, element_face)),
    element_face_unit_normal_->at(reflected_particle.element, element_face)
  );

  return reflected_particle;
}

std::shared_ptr<ParticleBoundary> ReflectingParticleBoundaryFactory::createBoundary(
  const Parameters& factory_params
) const {
  return std::make_shared<ReflectingParticleBoundary>(factory_params.element_face_unit_normal);
}

} // namespace mfpic
