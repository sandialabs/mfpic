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

  /// Checks whether the particles in the raw data are sorted first by element, and then by species.
  bool isSortedByElementThenSpecies() const;

  /**
   * @brief Get a subview of particles in a given element with a given species.
   *
   * Requires particles to be sorted with @ref sortByElementThenSpecies.
   * Since this is called in parallel, performance-critical kernels,
   * the status of the sort is not checked in optimized builds.
   *
   * @param[in] element Element in which to look for particles.
   * @param[in] species Species of particles for which to look.
   *
   * @returns A view of the requested particles, amenable to range-based for loops.
   */
  decltype(auto) particlesWithElementAndSpecies(this auto& self, int element, const Species& species) {
    assert(self.isSortedByElementThenSpecies());
    assert(element < (std::ssize(self.element_species_bins_) - 1) / self.numSpecies());

    const int species_index = self.getSpeciesIndex_(species);
    const int flattened_bin_index = element * self.numSpecies() + species_index;
    const int element_species_bin_begin = self.element_species_bins_[flattened_bin_index];
    const int element_species_bin_end = self.element_species_bins_[flattened_bin_index + 1];
    return std::span(
      std::next(self.particle_list_.begin(), element_species_bin_begin),
      std::next(self.particle_list_.begin(), element_species_bin_end)
    );
  }

private:
  /**
   * @brief Find the index of the given species in the internal species list.
   *
   * @param[in] species Species whose index to get.
   *
   * @returns Index of the given species in the internal species list.
   */
  int getSpeciesIndex_(const Species& species) const;

  /// List of particles.
  std::vector<Particle> particle_list_;

  /// List of species to which a particle has belonged at one point.
  std::vector<Species> species_represented_by_these_particles_;

  /// If container is sorted, indices representing bin edges in @a particle_list_ in which all particles have the same element and species.
  std::vector<int> element_species_bins_;
};

} // namespace mfpic
