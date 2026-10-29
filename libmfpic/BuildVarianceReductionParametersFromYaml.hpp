#pragma once

#include <string>
#include <vector>

namespace YAML {
class Node;
}

namespace mfpic {

/**
 * @brief Holds development options for variance reduction
 */
struct VarianceReductionParameters
{
  enum class Strategy {
    None,
    EulerFluid,
    UniformMaxwellian,
  };

  Strategy strategy = Strategy::None;
  /// Reference moments, used only by Strategy::UniformMaxwellian
  double reference_number_density = 0.0;
  std::vector<double> reference_bulk_velocity{0.0,0.0,0.0};
  double reference_temperature = 0.0;
  std::vector<std::string> reference_species;
  bool limit_variance_reduction = false;
  bool use_variance_reduced_electric_field = false;
};

/**
 * @brief Construct variance reduction parameters from yaml
 *
 * @param node - yaml node
 * @return VarianceReductionParameters
 */
VarianceReductionParameters buildVarianceReductionParametersFromYAML(const YAML::Node& node);

}
