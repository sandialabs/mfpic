#pragma once

#include <libmfpic/Species.hpp>

#include <mfem/mfem.hpp>
#include <yaml-cpp/yaml.h>

namespace mfpic {

/**
 * @brief parameters that define the state of a species at a single location.
 */
struct SourceStateParameters {
  double number_density;
  mfem::Vector bulk_velocity{0., 0., 0.};
  double temperature;
};

/**
 * @brief parameters that define a species over the whole domain, can either be used to load a species as a source or as initial
 *  conditions.
 */
struct SourceParameters {
  Species species;
  int num_particles = 0;

  SourceParameters() {};
  SourceParameters(
    const Species& species_in,
    const int num_particles_in = 0)
    : species(species_in)
    , num_particles(num_particles_in){};

  virtual ~SourceParameters() = default;

  /**
   * @brief Get an mfem::VectorCoefficient that represents an Euler fluid with the parameters in this object
   * 
   * @return std::unique_ptr<mfem::VectorCoefficient>
   */
  virtual std::unique_ptr<mfem::VectorCoefficient> getEulerVectorCoefficient() const = 0;
};

/**
 * @brief parameters defining a species that is constant in space
 */
struct ConstantSourceParameters : public SourceParameters {
  SourceStateParameters constant_state;

  ConstantSourceParameters(
    const Species& species,
    const SourceStateParameters& state_parameters,
    const int num_particles = 0);

  ConstantSourceParameters(
    const Species& species,
    const double number_density,
    const double temperature,
    const mfem::Vector& bulk_velocity=mfem::Vector({0., 0., 0.}));

  /**
   * @brief Get an mfem::VectorCoefficient that represents an constant euler fluid with parameters in this struct
   * 
   * @return std::unique_ptr<mfem::VectorCoefficient>
   */
  std::unique_ptr<mfem::VectorCoefficient> getEulerVectorCoefficient() const override;
};

/**
 * @brief parameters defining a species that is piecewise constant with one discontinuity 
 *  if x < discontinuity_location then species is given by left state otherwise species is right state
 */
struct SodSourceParameters : public SourceParameters {
  double discontinuity_location;
  SourceStateParameters left_state;
  SourceStateParameters right_state;

  SodSourceParameters(
    const Species& species,
    const double discontinuity_location,
    const SourceStateParameters& left_state_parameters,
    const SourceStateParameters& right_state_parameters,
    const int num_particles = 0);

  /**
   * @brief Get an mfem::VectorCoefficient that represents a piecewise constant euler fluid with parameters in this struct
   * 
   * @return std::unique_ptr<mfem::VectorCoefficient>
   */
  std::unique_ptr<mfem::VectorCoefficient> getEulerVectorCoefficient() const override;
};

/**
 * @brief build SourceStateParameters that defines a state from YAML
 * 
 * @param state_node - the YAML defining the state
 * @return SourceStateParameters - the parameters defining the state for a source or initial condition
 */
SourceStateParameters buildSourceStateParametersFromYAML(const YAML::Node& state_node);

/**
 * @brief build a list of parameters that define sources or initial conditions from YAML
 * 
 * @param sources_node - the YAML defining the sources or initial conditions
 * @param species_map - a map of species name to species definition
 * @return std::vector<std::unique_ptr<SourceParameters>> - a list of parameters defining the sources or initial conditions
 */
std::vector<std::unique_ptr<SourceParameters>> buildListOfSourceParametersFromYAML(
  const YAML::Node& sources_node,
  const std::unordered_map<std::string, Species>& species_map);

}