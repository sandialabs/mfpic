#include <libmfpic/ParticleBoundary.hpp>

namespace mfpic {

ParticleBoundary::~ParticleBoundary() = default;

void ParticleBoundaryFactory::setBoundaryAttribute(int boundary_attribute) {
  boundary_attribute_ = boundary_attribute;
}

int ParticleBoundaryFactory::getBoundaryAttribute() const {
  return boundary_attribute_;
}

ParticleBoundaryFactory::~ParticleBoundaryFactory() = default;

} // namespace mfpic
