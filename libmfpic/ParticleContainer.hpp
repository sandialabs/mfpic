#pragma once

#include <libmfpic/Particle.hpp>

namespace mfpic {

/// Enriched particle container.
class ParticleContainer {
public:
  using ParticleListType = std::vector<Particle>;

  /// Returns an iterator to the beginning of the particle container.
  ParticleListType::iterator begin();

  /// Returns an iterator to the beginning of the particle container.
  ParticleListType::const_iterator begin() const;

  /// Returns an iterator to the end of the particle container.
  ParticleListType::iterator end();

  /// Returns an iterator to the end of the particle container.
  ParticleListType::const_iterator end() const;

  /// Number of particles contained, both alive and dead.
  int numParticles() const;

  /// Number of species represented by the contained particles.
  int numSpecies() const;

  /**
   * @brief Add a single particle into the container.
   *
   * @param[in] particle Particle to add.
   */
  void addParticle(Particle particle);

  /**
   * @brief Add particles from another container into this one.
   *
   * @param[in] particles Particles to add.
   */
  void addParticles(const ParticleContainer& particles);

  /**
   * @brief Get the index of the given particle's species in this container's internal species list.
   *
   * @throws If the given particle's species is not in this container's internal species list.
   *
   * @param[in] particle Particle whose species' index to get.
   *
   * @returns The index of the given particle's species in this container's internal species list.
   */
  int getParticleSpeciesIndex(const Particle& particle) const;

  /// Remove all dead particles from the container.
  void cleanOutDeadParticles();

  /// Sort particles in the raw data first by element, and then by species.
  void sortByElementThenSpecies();

  /**
   * @brief Get a subview of particles in a given element with a given species.
   *
   * This internally sorts the particles by element and species so that the span is guaranteed contiguous,
   * but as a result it is not constant.
   *
   * @param[in] element Element in which to look for particles.
   * @param[in] species Species of particles for which to look.
   *
   * @returns A view of the requested particles, amenable to range-based for loops.
   */
  std::span<Particle> particlesWithElementAndSpecies(int element, const Species& species);

private:
  /// List of particles.
  std::vector<Particle> particle_list_;

  /// List of species to which a particle has belonged at one point.
  std::vector<Species> species_represented_by_these_particles_;

  /// If container is sorted, indices representing bin edges in @a particle_list_ in which all particles have the same element and species.
  std::vector<int> element_species_bins_;
};

} // namespace mfpic
