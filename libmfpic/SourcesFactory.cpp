#include <libmfpic/Errors.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/SourcesFactory.hpp>

namespace mfpic {

namespace {
mfem::Vector constructPrimitiveState(const SourceStateParameters state) {
  return euler::constructPrimitiveState(state.number_density, state.bulk_velocity, state.temperature);
}
}

ConstantSourceParameters::ConstantSourceParameters(
  const Species& species,
  const SourceStateParameters& state_parameters,
  const int num_particles)
  : SourceParameters(species, num_particles)
  , constant_state(state_parameters)
  {}

ConstantSourceParameters::ConstantSourceParameters(
  const Species& species,
  const double number_density,
  const double temperature,
  const mfem::Vector& bulk_velocity)
  : SourceParameters(species)
  , constant_state{.number_density = number_density, .bulk_velocity = bulk_velocity, .temperature = temperature}
{}

std::unique_ptr<mfem::VectorCoefficient> ConstantSourceParameters::getEulerVectorCoefficient() const {
  const mfem::Vector primitive_state = constructPrimitiveState(constant_state);

  const mfem::Vector conservative_state = euler::convertFromPrimitiveToConservative(primitive_state, species);
  auto euler_coefficient = std::make_unique<mfem::VectorConstantCoefficient>(conservative_state);
  return euler_coefficient;
}

SodSourceParameters::SodSourceParameters(
  const Species& species,
  const double discontinuity_location_in,
  const SourceStateParameters& left_state_parameters,
  const SourceStateParameters& right_state_parameters,
  const int num_particles)
  : SourceParameters(species, num_particles)
  , discontinuity_location(discontinuity_location_in)
  , left_state(left_state_parameters)
  , right_state(right_state_parameters)
{}

std::unique_ptr<mfem::VectorCoefficient> SodSourceParameters::getEulerVectorCoefficient() const {
  const mfem::Vector primitive_state_left = constructPrimitiveState(left_state);
  const mfem::Vector primitive_state_right = constructPrimitiveState(right_state);

  const mfem::Vector conservative_state_left = euler::convertFromPrimitiveToConservative(primitive_state_left, species);
  const mfem::Vector conservative_state_right = euler::convertFromPrimitiveToConservative(primitive_state_right, species);

  auto function = [
    &discontinuity_location = discontinuity_location,
    conservative_state_right = conservative_state_right,
    conservative_state_left = conservative_state_left](const mfem::Vector &x, mfem::Vector &y)
  {
    if (x[0] < discontinuity_location) {
      y = conservative_state_left;
    } else {
      y = conservative_state_right;
    }
  };
  auto euler_coefficient = std::make_unique<mfem::VectorFunctionCoefficient>(euler::ConservativeVariables::NUM_VARS, function);
  return euler_coefficient;
}

SourceStateParameters buildSourceStateParametersFromYAML(const YAML::Node& state_node) {
  const double number_density = state_node["Number Density"].as<double>();
  if (number_density <= 0.0) {
    errorWithUserMessage(formatParseMessage(state_node["Number Density"], "Number Density is nonpositive!"));
  }

  const YAML::Node& bulk_velocity_node = state_node["Bulk Velocity"];
  mfem::Vector bulk_velocity({0.0, 0.0, 0.0});
  if (bulk_velocity_node) {
    if (bulk_velocity_node.IsSequence()) {
      if (bulk_velocity_node.size() != 3) {
        errorWithUserMessage(formatParseMessage(bulk_velocity_node, "Bulk Velocity does not have length 3!"));
      }
      for (int i = 0; i < 3; i++) {
        bulk_velocity[i] = bulk_velocity_node[i].as<double>();
      }
    } else {
      std::cout << formatParseMessage(
        bulk_velocity_node, "Bulk Velocity is not a sequence, so it's ignored. Continuing.") << std::endl;
    }
  }

  const double temperature = state_node["Temperature"].as<double>();
  if (temperature < 0.0) {
    errorWithUserMessage(formatParseMessage(state_node["Temperature"], "Temperature is negative!"));
  }

  SourceStateParameters state_parameters{
    .number_density = number_density,
    .bulk_velocity = bulk_velocity,
    .temperature = temperature};
  return state_parameters;
}

std::vector<std::unique_ptr<SourceParameters>> buildListOfSourceParametersFromYAML(
  const YAML::Node& sources_node,
  const std::unordered_map<std::string, Species>& species_map)
{
  assert(sources_node.IsSequence() or sources_node.IsNull() or not sources_node.IsDefined());

  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters;

  if (sources_node.IsSequence()) {
    for (const YAML::Node& source : sources_node) {
      std::vector<std::string> species_names;
      const YAML::Node& species_node = source["Species"];
      if (species_node.IsSequence()) {
        for (const YAML::Node& species : species_node) {
          const std::string species_name = species.as<std::string>();
          if (not species_map.contains(species_name)) {
            errorWithUserMessage(formatParseMessage(species, "Species was not created in the Species block!"));
          }
          species_names.push_back(species_name);
        }
      } else {
        errorWithUserMessage(formatParseMessage(species_node, "Species is not a sequence!"));
      }

      int num_particles_per_species = 0;
      if (source["Number of Macroparticles per Species"]) {
        num_particles_per_species = source["Number of Macroparticles per Species"].as<int>();
        if (num_particles_per_species <= 0) {
          errorWithUserMessage(formatParseMessage(
            source["Number of Macroparticles per Species"], "Number of Macroparticles per Species is nonpositive!"));
        }
      }

      if (source["Constant"]) {
        const YAML::Node& state_node = source["Constant"];
        const SourceStateParameters state_parameters = buildSourceStateParametersFromYAML(state_node);
        for (const std::string& species_name : species_names) {
          list_of_parameters.push_back(
            std::make_unique<ConstantSourceParameters>(
              species_map.at(species_name),
              state_parameters,
              num_particles_per_species));
        }
      } else if (source["Sod"]) {
        const YAML::Node& sod_node = source["Sod"];
        const double discontinuity_location = sod_node["Discontinuity Location"].as<double>();
        const YAML::Node& left_state_node = sod_node["Left State"];
        const SourceStateParameters left_state_parameters = buildSourceStateParametersFromYAML(left_state_node);

        const YAML::Node& right_state_node = sod_node["Right State"];
        const SourceStateParameters right_state_parameters = buildSourceStateParametersFromYAML(right_state_node);

        for (const std::string& species_name : species_names) {
          list_of_parameters.push_back(std::make_unique<SodSourceParameters>(species_map.at(species_name), discontinuity_location, left_state_parameters, right_state_parameters, num_particles_per_species));
        }

      } else {
        errorWithUserMessage(formatParseMessage(source, "It is required to either specify \"Constant\" or \"Sod\"."));
      }
    }
  } else {
    std::cout << formatParseMessage(
      sources_node, "Not a sequence, so it's ignored. Continuing.") << std::endl;
  }

  return list_of_parameters;
}

}