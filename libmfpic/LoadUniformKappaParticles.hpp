#pragma once

#include <libmfpic/GenerateKappaVelocity.hpp>
#include <libmfpic/MeshUtilities.hpp>
#include <libmfpic/ParticleContainer.hpp>
#include <libmfpic/Species.hpp>
#include <libmfpic/UniformMeshDistribution.hpp>

#include <mfem/mfem.hpp>

#include <random>

namespace mfpic {

/**
 * @brief Load a requested number of particles distributed uniformly in space and Kappa distribution in velocity space.
 *
 * @tparam Generator A UniformRandomBitGenerator type.
 *
 * @param[in]     species_list               List of species for which to create particles.
 * @param[in]     bulk_velocity              Bulk velocity of particle distribution.
 * @param[in]     temperature                Species-wise temperature of particle distribution.
 * @param[in]     kappa                      Kappa parameter in the kappa-distribution
 * @param[in]     number_density_per_species Number density of each species.
 * @param[in]     num_particles_per_species  Number of particles to create for each species.
 * @param[in,out] generator                  A UniformRandomBitGenerator used to generate some random numbers.
 * @param[in]     mesh                       Mesh in which to create particles.
 *
 * @returns Container of created particles.
 */
template <std::uniform_random_bit_generator Generator>
ParticleContainer loadUniformKappaParticles(
  Species species,
  mfem::Vector bulk_velocity,
  double temperature,
  double kappa,
  double number_density,
  int num_particles,
  Generator& generator,
  std::shared_ptr<mfem::Mesh> mesh
) {
  ParticleContainer particles;

  const double mesh_volume = getMeshVolume(*mesh);
  const double particle_weight = number_density * mesh_volume / num_particles;
  UniformMeshDistribution position_distribution(mesh);

  for (int i = 0; i < num_particles; ++i) {
    const mfem::Vector velocity = generateKappaVelocity(bulk_velocity, temperature, kappa, species.mass, generator);

    mfem::Vector position({0.0, 0.0, 0.0});
    const auto [random_mesh_position, element] = position_distribution.generateRandomPointAndElement(generator);
    for (int i = 0; i < random_mesh_position.Size(); i++) {
      position[i] = random_mesh_position[i];
    }
    particles.addParticle(Particle{
      .position = position,
      .velocity = velocity,
      .element = element,
      .species = species,
      .weight = particle_weight,
      .is_alive = true,
    });
  }

  return particles;
}

} // namespace mfpic
