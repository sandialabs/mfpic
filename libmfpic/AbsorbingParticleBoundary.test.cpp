#include <libmfpic/AbsorbingParticleBoundary.hpp>
#include <libmfpic/Particle.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

TEST(AbsorbingParticleBoundary, AbsorbedParticleDies) {
  AbsorbingParticleBoundary boundary;
  Particle incoming_particle;

  constexpr int unused_face = 0;
  Particle absorbed_particle = boundary.applyBoundary(unused_face, incoming_particle);

  ASSERT_FALSE(absorbed_particle.is_alive);
}

} // namespace
