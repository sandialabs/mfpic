#include <libmfpic/ParticleBoundary.hpp>

namespace mfpic {

void ParticleBoundaryFactory::setBoundaryAttribute(int boundary_attribute) {
  boundary_attribute_ = boundary_attribute;
}

int ParticleBoundaryFactory::getBoundaryAttribute() const {
  return boundary_attribute_;
}

} // namespace mfpic
