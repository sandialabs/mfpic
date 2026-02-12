#include <libmfpic/ParticleContainer.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

TEST(ParticleContainer, DefaultInitializedParticleContainerIsEmpty) {
  ParticleContainer particle_container;

  EXPECT_EQ(particle_container.numParticles(), 0);
}

TEST(ParticleContainer, AddParticleAddsAParticle) {
  ParticleContainer particle_container;
  constexpr int element = 20;
  Particle particle_to_add{
    .element = element,
  };
  particle_container.addParticle(particle_to_add);

  EXPECT_EQ(particle_container.numParticles(), 1);
  const Particle& added_particle = *particle_container.begin();
  EXPECT_EQ(added_particle.element, element);
}

TEST(ParticleContainer, AddParticlesAddsParticles) {
  ParticleContainer particle_container_to_add_into;
  ParticleContainer particle_container_to_add_from;
  constexpr int element = 20;
  Particle particle_to_add{
    .element = element,
  };
  particle_container_to_add_from.addParticle(particle_to_add);

  particle_container_to_add_into.addParticles(particle_container_to_add_from);

  EXPECT_EQ(particle_container_to_add_into.numParticles(), 1);
  const Particle& added_particle = *particle_container_to_add_into.begin();
  EXPECT_EQ(added_particle.element, element);
}

} // namespace
