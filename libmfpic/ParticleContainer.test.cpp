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

void checkThatParticlesAreSortedByElementThenSpecies(const ParticleContainer& particles) {
  for (
    auto particle_iter = particles.begin();
    std::distance(particle_iter, particles.end()) > 1;
    std::advance(particle_iter, 1)
  ) {
    const Particle& this_particle = *particle_iter;
    const Particle& next_particle = *std::next(particle_iter);
    if (this_particle.element == next_particle.element) {
      ASSERT_LE(
        particles.getParticleSpeciesIndex(this_particle),
        particles.getParticleSpeciesIndex(next_particle)
      );
    }
    else {
      ASSERT_LT(this_particle.element, next_particle.element);
    }
  }
}

TEST(ParticleContainer, SortingAlreadySortedParticlesPreservesSortedness) {
  ParticleContainer already_sorted;
  already_sorted.addParticle(Particle{.element = 0, .species = species_0});
  already_sorted.addParticle(Particle{.element = 0, .species = species_1});
  already_sorted.addParticle(Particle{.element = 1, .species = species_0});
  already_sorted.addParticle(Particle{.element = 1, .species = species_1});

  already_sorted.sortByElementThenSpecies();

  checkThatParticlesAreSortedByElementThenSpecies(already_sorted);
}

TEST(ParticleContainer, SortByElementsWithOneSpecies) {
  ParticleContainer to_sort;
  to_sort.addParticle(Particle{.element = 3});
  to_sort.addParticle(Particle{.element = 0});
  to_sort.addParticle(Particle{.element = 1});

  to_sort.sortByElementThenSpecies();

  checkThatParticlesAreSortedByElementThenSpecies(to_sort);
}

TEST(ParticleContainer, SortBySpeciesWithAllParticlesInElementZero) {
  ParticleContainer to_sort;
  to_sort.addParticle(Particle{.element = 0, .species = species_2});
  to_sort.addParticle(Particle{.element = 0, .species = species_1});
  to_sort.addParticle(Particle{.element = 0, .species = species_0});

  to_sort.sortByElementThenSpecies();

  checkThatParticlesAreSortedByElementThenSpecies(to_sort);
}

TEST(ParticleContainer, SortBySpeciesWithAllParticlesInElementNonzero) {
  constexpr int element = 3;
  ParticleContainer to_sort;
  to_sort.addParticle(Particle{.element = element, .species = species_2});
  to_sort.addParticle(Particle{.element = element, .species = species_1});
  to_sort.addParticle(Particle{.element = element, .species = species_0});

  to_sort.sortByElementThenSpecies();

  checkThatParticlesAreSortedByElementThenSpecies(to_sort);
}

TEST(ParticleContainer, SortByElementThenSpecies) {
  ParticleContainer to_sort;
  to_sort.addParticle(Particle{.element = 4, .species = species_1});
  to_sort.addParticle(Particle{.element = 2, .species = species_1});
  to_sort.addParticle(Particle{.element = 2, .species = species_0});
  to_sort.addParticle(Particle{.element = 2, .species = species_0});
  to_sort.addParticle(Particle{.element = 3, .species = species_0});

  to_sort.sortByElementThenSpecies();

  checkThatParticlesAreSortedByElementThenSpecies(to_sort);
}

} // namespace
