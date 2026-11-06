#include <libmfpic/BuildVarianceReductionParametersFromYaml.hpp>
#include <libmfpic/Errors.hpp>
#include <libmfpic/ParseHelpers.hpp>

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
    const mfem::Vector bulk_velocity = parseVector3(moments_node["Bulk Velocity"], "Bulk Velocity");
    for (int i = 0; i < 3; ++i) {
      parameters.reference_bulk_velocity[i] = bulk_velocity[i];
    }
    parameters.reference_temperature = moments_node["Temperature"].as<double>();
    if (moments_node["Species"]) {
      parameters.reference_species = moments_node["Species"].as<std::vector<std::string>>();
    }
  } else if (moments_node) {
    errorWithUserMessage(formatParseMessage(moments_node,
      "Uniform Maxwellian Moments is only used with Strategy: Uniform Maxwellian."));
  }

  const YAML::Node& cell_limiting_node = node["Cell Limiting"];
  if (cell_limiting_node) {
    const std::string cell_limiting = cell_limiting_node.as<std::string>();
    if (cell_limiting == "None") {
      parameters.cell_limiting = VarianceReductionParameters::CellLimiting::None;
    } else if (cell_limiting == "Every Particle") {
      parameters.cell_limiting = VarianceReductionParameters::CellLimiting::EveryParticle;
    } else if (cell_limiting == "Variance") {
      parameters.cell_limiting = VarianceReductionParameters::CellLimiting::Variance;
    } else {
      errorWithUserMessage(formatParseMessage(cell_limiting_node,
        "Unknown Cell Limiting \"" + cell_limiting + "\". Must be \"None\", \"Every Particle\", or \"Variance\"."));
    }
  }

  parameters.use_variance_reduced_electric_field = node["Use VR E Field"].as<bool>(false);

  return parameters;
}

} // namespace mfpic
