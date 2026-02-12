#include <libmfpic/AbsorbingParticleBoundary.hpp>
#include <libmfpic/Particle.hpp>

namespace mfpic {

AbsorbingParticleBoundary::AbsorbingParticleBoundary() = default;

Particle AbsorbingParticleBoundary::applyBoundary(int, const Particle &incoming_particle) const {
  Particle dead_particle = incoming_particle;
  dead_particle.is_alive = false;
  return dead_particle;
}

std::shared_ptr<ParticleBoundary> AbsorbingParticleBoundaryFactory::createBoundary(const Parameters&) const {
  return std::make_shared<AbsorbingParticleBoundary>();
}

} // namespace mfpic
