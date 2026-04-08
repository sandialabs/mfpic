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
  const SourceParameters& source_parameters,
  Generator& generator,
  std::shared_ptr<mfem::Mesh> mesh
) {
  ParticleContainer particles;

  auto number_density_function = [&] (const mfem::Vector& x) -> double {
    SourceStateParameters params_at_point = source_parameters.sourceStateParametersAtPoint(x);
    return params_at_point.number_density;
  };
  const mfem::Vector physical_particles_per_element = elementwiseIntegral(*mesh, number_density_function);
  const double total_physical_particles = std::reduce(
    physical_particles_per_element.begin(),
    physical_particles_per_element.end()
  );
  const double particle_weight = total_physical_particles / source_parameters.num_particles;

  MeshDistribution position_distribution(mesh, number_density_function);

  const Species& species = source_parameters.species;
  for (int i = 0; i < source_parameters.num_particles; ++i) {
    mfem::Vector position({0.0, 0.0, 0.0});
    const auto [random_mesh_position, element] = position_distribution.generateRandomPointAndElement(generator);
    for (int i = 0; i < random_mesh_position.Size(); i++) {
      position[i] = random_mesh_position[i];
    }

    SourceStateParameters source_state_parameters = source_parameters.sourceStateParametersAtPoint(position);

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
