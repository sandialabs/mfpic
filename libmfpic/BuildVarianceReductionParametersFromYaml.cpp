#include <libmfpic/BuildVarianceReductionParametersFromYaml.hpp>
#include <libmfpic/Errors.hpp>

#include <yaml-cpp/yaml.h>

namespace mfpic {

VarianceReductionParameters buildVarianceReductionParametersFromYAML(const YAML::Node& node)
{
  VarianceReductionParameters parameters;

  const std::string strategy = node["Strategy"].as<std::string>("Euler Fluid");
  if (strategy == "none") parameters.strategy = VarianceReductionParameters::Strategy::None;
  else if (strategy == "Euler Fluid") parameters.strategy = VarianceReductionParameters::Strategy::EulerFluid;
  else if (strategy == "Local Maxwellian") parameters.strategy = VarianceReductionParameters::Strategy::LocalMaxwellian;
  else
    throw std::runtime_error("Unknown variance reduction Strategy: " + strategy);

  if (node["Specified LF Moments"]) {
    parameters.specified_lf_moments = true;
    const YAML::Node& specified_moments_node = node["Specified LF Moments"];
    parameters.reference_number_density = specified_moments_node["Reference Number Density"].as<double>();
    const auto reference_bulk_velocity = specified_moments_node["Reference Bulk Velocity"].as<std::vector<double>>();
    for (std::size_t i = 0; i < std::min<std::size_t>(3, reference_bulk_velocity.size()); ++i)
    {
      parameters.reference_bulk_velocity[i] = reference_bulk_velocity[i];
    }
    parameters.reference_temperature = specified_moments_node["Reference Temperature"].as<double>();
  }

  parameters.limit_variance_reduction = node["Apply Cell Limiting"].as<bool>(true);

  parameters.use_variance_reduced_electric_field = node["Use VR E Field"].as<bool>(false);

  return parameters;
} 

} // namespace mfpic
