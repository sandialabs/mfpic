#include <libmfpic/BuildVarianceReductionParametersFromYaml.hpp>
#include <libmfpic/Errors.hpp>

#include <yaml-cpp/yaml.h>

namespace mfpic {

VarianceReductionParameters buildVarianceReductionParametersFromYAML(const YAML::Node& node)
{
  VarianceReductionParameters parameters;

  const YAML::Node& strategy_node = node["Strategy"];
  const std::string strategy = strategy_node.as<std::string>("Euler Fluid");
  if (strategy == "none") parameters.strategy = VarianceReductionParameters::Strategy::None;
  else if (strategy == "Euler Fluid") parameters.strategy = VarianceReductionParameters::Strategy::EulerFluid;
  else if (strategy == "Uniform Maxwellian") parameters.strategy = VarianceReductionParameters::Strategy::UniformMaxwellian;
  else {
    errorWithUserMessage(formatParseMessage(strategy_node,
      "Unknown variance reduction Strategy \"" + strategy + "\". Must be \"none\", \"Euler Fluid\", or \"Uniform Maxwellian\"."));
  }

  const YAML::Node& moments_node = node["Uniform Maxwellian Moments"];
  if (parameters.strategy == VarianceReductionParameters::Strategy::UniformMaxwellian) {
    if (not moments_node) {
      errorWithUserMessage(formatParseMessage(node, "Strategy: Uniform Maxwellian requires a Uniform Maxwellian Moments block."));
    }
    parameters.reference_number_density = moments_node["Number Density"].as<double>();
    const YAML::Node& bulk_velocity_node = moments_node["Bulk Velocity"];
    if (not bulk_velocity_node.IsSequence() or bulk_velocity_node.size() != 3) {
      errorWithUserMessage(formatParseMessage(bulk_velocity_node, "Bulk Velocity does not have length 3!"));
    }
    for (int i = 0; i < 3; ++i) {
      parameters.reference_bulk_velocity[i] = bulk_velocity_node[i].as<double>();
    }
    parameters.reference_temperature = moments_node["Temperature"].as<double>();
    if (moments_node["Species"]) {
      parameters.reference_species = moments_node["Species"].as<std::vector<std::string>>();
    }
  } else if (moments_node) {
    errorWithUserMessage(formatParseMessage(moments_node,
      "Uniform Maxwellian Moments is only used with Strategy: Uniform Maxwellian."));
  }
  if (node["Specified LF Moments"]) {
    errorWithUserMessage(formatParseMessage(node["Specified LF Moments"],
      "Specified LF Moments has been replaced by Strategy: Uniform Maxwellian with a Uniform Maxwellian Moments block."));
  }

  parameters.limit_variance_reduction = node["Apply Cell Limiting"].as<bool>(false);

  parameters.use_variance_reduced_electric_field = node["Use VR E Field"].as<bool>(false);

  return parameters;
}

} // namespace mfpic
