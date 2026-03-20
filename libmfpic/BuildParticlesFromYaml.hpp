#pragma once

#include <libmfpic/Errors.hpp>
#include <libmfpic/LoadUniformMaxwellianParticles.hpp>
#include <libmfpic/LoadUniformKappaParticles.hpp>
#include <libmfpic/ParticleContainer.hpp>
#include <libmfpic/SourcesFactory.hpp>

#include <yaml-cpp/yaml.h>

namespace mfpic {

/**
 * @brief Construct particles in the given mesh using a specification from a YAML file.
 *
 * @tparam Generator A UniformRandomBitGenerator type.
 *
 * @param[in]     sources_node   A YAML node containing a sequence of source parameters.
 * @param[in]     species_map    Map species names to Species structs.
 * @param[in,out] generator      A UniformRandomBitGenerator used to generate some random numbers.
 * @param[in]     mesh           Mesh in which to generate particles.
 */
template <std::uniform_random_bit_generator Generator>
ParticleContainer buildParticlesFromYaml(
  const YAML::Node& sources_node,
  const std::unordered_map<std::string, Species>& species_map,
  Generator& generator,
  std::shared_ptr<mfem::Mesh> mesh
);

// template definitions

template <std::uniform_random_bit_generator Generator>
ParticleContainer buildParticlesFromYaml(
  const YAML::Node& sources_node,
  const std::unordered_map<std::string, Species>& species_map,
  Generator& generator,
  std::shared_ptr<mfem::Mesh> mesh
) {
  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters = buildListOfSourceParametersFromYAML(
    sources_node,
    species_map
  );

  ParticleContainer particles;
  for (const std::unique_ptr<SourceParameters>& parameters : list_of_parameters) {
    auto constant_parameters = dynamic_cast<const ConstantSourceParameters&>(*parameters);

    if (constant_parameters.constant_state.kappa == -1)
    {
      particles.addParticles(loadUniformMaxwellianParticles(
        parameters->species,
        constant_parameters.constant_state.bulk_velocity,
        constant_parameters.constant_state.temperature,
        constant_parameters.constant_state.number_density,
        parameters->num_particles,
        generator,
        mesh
      ));
    }
    else
    {
      particles.addParticles(loadUniformKappaParticles(
        parameters->species,
        constant_parameters.constant_state.bulk_velocity,
        constant_parameters.constant_state.temperature,
        constant_parameters.constant_state.kappa,
        constant_parameters.constant_state.number_density,
        parameters->num_particles,
        generator,
        mesh
      ));
    }
  }

  return particles;
}

} // namespace mfpic
