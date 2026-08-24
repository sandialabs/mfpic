#pragma once

#include <libmfpic/Errors.hpp>
#include <libmfpic/LoadParticles.hpp>
#include <libmfpic/ParticleContainer.hpp>
#include <libmfpic/SourcesFactory.hpp>

#include <yaml-cpp/yaml.h>

namespace mfpic {

double evaluateSingleParticleDistributionFunction(
  const SourceParameters& source_parameters,
  const mfem::Vector& x,
  const mfem::Vector& v
) {
  const SourceStateParameters ssp = source_parameters.sourceStateParametersAtPoint(x);
  const Species& species = source_parameters.species;

  mfem::Vector primitive_state(5);
  primitive_state(0) = ssp.number_density;
  primitive_state(1) = ssp.bulk_velocity(0);
  primitive_state(2) = ssp.bulk_velocity(1);
  primitive_state(3) = ssp.bulk_velocity(2);
  primitive_state(4) = ssp.temperature;

  if (ssp.kappa > 0.0) {
    return euler::evaluateIsotropicKappaDistribution(primitive_state, v, ssp.kappa, species,3);
  } else {
    return euler::evaluateMaxwellian(primitive_state, v, species);
  }
}

void updateParticleDistributionFunctionValues(
  ParticleContainer& particles,
  const std::vector<std::unique_ptr<SourceParameters>>&list_of_parameters
) {
  for (auto& particle : particles) { 
    double f_global = 0.0;

    for (const auto& source_parameters: list_of_parameters) {

      if (!(particle.species==source_parameters->species)) {
        continue;
      }
      f_global += evaluateSingleParticleDistributionFunction(*source_parameters, particle.position, particle.velocity);
    }

    particle.particle_distribution_function_value = f_global;
  }
}

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
  std::shared_ptr<mfem::Mesh> mesh,
  int velocity_dims=3
) {
  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters = buildListOfSourceParametersFromYAML(
    sources_node,
    species_map
  );

  ParticleContainer particles;
  for (const std::unique_ptr<SourceParameters>& parameters : list_of_parameters) {
    particles.addParticles(loadParticles(
      *parameters,
      generator,
      mesh,
      velocity_dims
    ));
  }

  //updateParticleDistributionFunctionValues(particles, list_of_parameters);

  return particles;
}


} // namespace mfpic
