#pragma once

#include <libmfpic/Constants.hpp>
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
  assert(bulk_velocity.Size() == 3);
  assert(temperature >= 0.0);

  ParticleContainer particles;

  const double mesh_volume = getMeshVolume(*mesh);
  const double particle_weight = number_density * mesh_volume / num_particles;
  UniformMeshDistribution position_distribution(mesh);

  //TODO: Take a second look at 2 here versus reference
  const double v_thermal_squared = 2.0 * constants::boltzmann_constant * temperature / species.mass;
  for (int i = 0; i < num_particles; ++i) {
    mfem::Vector velocity = bulk_velocity;
    if (v_thermal_squared > 0.0) {
      const double nu = 2.0 * kappa + 1.0;
      const double scale = std::sqrt(kappa * v_thermal_squared / (2.0 * kappa + 1.0));
      std::student_t_distribution<double> velocity_distribution(nu);
      for (int dimension = 0; dimension < 3; dimension++) {
        velocity[dimension] += scale * velocity_distribution(generator); 
      }
    }

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
