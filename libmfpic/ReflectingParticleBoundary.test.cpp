#include <libmfpic/ElementFaceContainer.hpp>
#include <libmfpic/Particle.hpp>
#include <libmfpic/ReflectingParticleBoundary.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

TEST(ReflectingParticleBoundary, HeadOnParticleReflectsStraightBack) {
  Particle particle{
    .velocity = mfem::Vector{1.0, 0.0, 0.0},
    .element = 0,
  };
  ElementFaceContainer<mfem::Vector> element_face_unit_normal;
  element_face_unit_normal.insert(0, 0, mfem::Vector{1.0, 0.0, 0.0});
  ReflectingParticleBoundary boundary(std::make_shared<ElementFaceContainer<mfem::Vector>>(element_face_unit_normal));

  particle = boundary.applyBoundary(0, particle);

  EXPECT_DOUBLE_EQ(particle.velocity[0], -1.0);
  EXPECT_DOUBLE_EQ(particle.velocity[1], 0.0);
  EXPECT_DOUBLE_EQ(particle.velocity[2], 0.0);
}

TEST(ReflectingParticleBoundary, ObliqueReflectionDoesNotAffectTangentVelocityComponents) {
  Particle particle{
    .velocity = mfem::Vector{-2342354.0, 5051231.0, 132213643.0},
    .element = 0,
  };
  ElementFaceContainer<mfem::Vector> element_face_unit_normal;
  element_face_unit_normal.insert(0, 0, mfem::Vector{-1.0, 0.0, 0.0});
  ReflectingParticleBoundary boundary(std::make_shared<ElementFaceContainer<mfem::Vector>>(element_face_unit_normal));

  particle = boundary.applyBoundary(0, particle);

  EXPECT_DOUBLE_EQ(particle.velocity[0], 2342354.0);
  EXPECT_DOUBLE_EQ(particle.velocity[1], 5051231.0);
  EXPECT_DOUBLE_EQ(particle.velocity[2], 132213643.0);
}

} // namespace
