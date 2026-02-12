#pragma once

#include <libmfpic/ParticleBoundary.hpp>

namespace mfpic {

/// Absorbs particles incident on faces.
class AbsorbingParticleBoundary : public ParticleBoundary {
public:
  AbsorbingParticleBoundary();

  virtual Particle applyBoundary(int, const Particle &incoming_particle) const override;
};

class AbsorbingParticleBoundaryFactory : public ParticleBoundaryFactory {
public:
  virtual std::shared_ptr<ParticleBoundary> createBoundary(
    const Parameters&
  ) const override;
};

} // namespace mfpic
