#include <libmfpic/ParticleContainer.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

const Species species_0{.name = "species_0"};
const Species species_1{.name = "species_1"};
const Species species_2{.name = "species_2"};

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

TEST(ParticleContainer, SortingAlreadySortedParticlesPreservesSortedness) {
  ParticleContainer already_sorted;
  already_sorted.addParticle(Particle{.element = 0, .species = species_0});
  already_sorted.addParticle(Particle{.element = 0, .species = species_1});
  already_sorted.addParticle(Particle{.element = 1, .species = species_0});
  already_sorted.addParticle(Particle{.element = 1, .species = species_1});

  ASSERT_TRUE(already_sorted.isSortedByElementThenSpecies());
  already_sorted.sortByElementThenSpecies();
  ASSERT_TRUE(already_sorted.isSortedByElementThenSpecies());
}

TEST(ParticleContainer, SortByElementsWithOneSpecies) {
  ParticleContainer to_sort;
  to_sort.addParticle(Particle{.element = 3});
  to_sort.addParticle(Particle{.element = 0});
  to_sort.addParticle(Particle{.element = 1});

  ASSERT_FALSE(to_sort.isSortedByElementThenSpecies());
  to_sort.sortByElementThenSpecies();
  ASSERT_TRUE(to_sort.isSortedByElementThenSpecies());
}

TEST(ParticleContainer, SortByElementThenSpecies) {
  ParticleContainer to_sort;
  to_sort.addParticle(Particle{.element = 3, .species = species_0});
  to_sort.addParticle(Particle{.element = 2, .species = species_0});
  to_sort.addParticle(Particle{.element = 3, .species = species_0});
  to_sort.addParticle(Particle{.element = 4, .species = species_1});
  to_sort.addParticle(Particle{.element = 2, .species = species_1});

  ASSERT_FALSE(to_sort.isSortedByElementThenSpecies());
  to_sort.sortByElementThenSpecies();
  ASSERT_TRUE(to_sort.isSortedByElementThenSpecies());
}

TEST(ParticleContainer, ParticlesWithElementAndSpeciesEmptyIfNoSuchParticlesExist) {
  ParticleContainer no_particles_in_element_1_with_species_0;
  no_particles_in_element_1_with_species_0.addParticle(Particle{.element = 0, .species = species_0});
  no_particles_in_element_1_with_species_0.addParticle(Particle{.element = 1, .species = species_1});

  no_particles_in_element_1_with_species_0.sortByElementThenSpecies();
  std::span<Particle> particles_in_element_1_with_species_0 =
    no_particles_in_element_1_with_species_0.particlesWithElementAndSpecies(1, species_0);

  ASSERT_EQ(0, particles_in_element_1_with_species_0.size());
}

TEST(ParticleContainer, SizeOfParticlesWithElementAndSpeciesMatchesExpectedNumber) {
  ParticleContainer particles;
  constexpr int expected_num_particles = 4;
  constexpr int checked_element = 3;
  const Species checked_species = species_0;
  for (int i = 0; i < expected_num_particles; i++) {
    particles.addParticle(Particle{.element = checked_element, .species = checked_species});
    constexpr int unchecked_element = 2;
    const Species unchecked_species = species_1;
    particles.addParticle(Particle{.element = unchecked_element, .species = unchecked_species});
  }

  particles.sortByElementThenSpecies();
  std::span<Particle> checked_particles = particles.particlesWithElementAndSpecies(checked_element, checked_species);

  ASSERT_EQ(expected_num_particles, checked_particles.size());
}

TEST(ParticleContainer, AllParticlesWithElementAndSpeciesActuallyHaveRequestedElementAndSpecies) {
  ParticleContainer particles;
  constexpr int expected_num_particles = 4;
  constexpr int requested_element = 3;
  const Species requested_species = species_0;
  for (int i = 0; i < expected_num_particles; i++) {
    particles.addParticle(Particle{.element = requested_element, .species = requested_species});
    constexpr int unrequested_element = 2;
    const Species unrequested_species = species_1;
    particles.addParticle(Particle{.element = unrequested_element, .species = unrequested_species});
  }

  particles.sortByElementThenSpecies();
  for (const Particle& particle : particles.particlesWithElementAndSpecies(requested_element, requested_species)) {
    EXPECT_EQ(particle.element, requested_element);
    EXPECT_EQ(particle.species, requested_species);
  }
}

} // namespace
