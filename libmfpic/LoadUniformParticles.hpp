#pragma once

#include <libmfpic/GenerateKappaVelocity.hpp>
#include <libmfpic/GenerateMaxwellianVelocity.hpp>
#include <libmfpic/MeshDistribution.hpp>
#include <libmfpic/MeshUtilities.hpp>
#include <libmfpic/ParticleContainer.hpp>
#include <libmfpic/SourcesFactory.hpp>
#include <libmfpic/Species.hpp>

#include <mfem/mfem.hpp>

#include <random>

namespace mfpic {

/**
 * @brief Load a requested number of particles distributed uniformly in space and a parametrized velocity distribuion.
 *
 * @tparam Generator A UniformRandomBitGenerator type.
 *
 * @param[in]     species                    Species for which to create particles.
 * @param[in]     source_state_parameters    Parameters for the state of the particle distribution.
 * @param[in]     num_particles              Number of particles to create.
 * @param[in,out] generator                  A UniformRandomBitGenerator used to generate some random numbers.
 * @param[in]     mesh                       Mesh in which to create particles.
 *
 * @returns Container of created particles.
 */
template <std::uniform_random_bit_generator Generator>
ParticleContainer loadUniformParticles(
  Species species,
  const SourceStateParameters& source_state_parameters,
  int num_particles,
  Generator& generator,
  std::shared_ptr<mfem::Mesh> mesh
) {
  ParticleContainer particles;

  const double mesh_volume = getMeshVolume(*mesh);
  const double particle_weight = source_state_parameters.number_density * mesh_volume / num_particles;
  MeshDistribution position_distribution(mesh);

  for (int i = 0; i < num_particles; ++i) {
    const double kappa = source_state_parameters.kappa;
    mfem::Vector velocity;
    if (kappa > 0.0) {
      velocity = generateKappaVelocity(
        source_state_parameters.bulk_velocity,
        source_state_parameters.temperature,
        kappa,
        species.mass,
        generator
      );
    }
    else {
      velocity = generateMaxwellianVelocity(
        source_state_parameters.bulk_velocity,
        source_state_parameters.temperature,
        species.mass,
        generator
      );
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
