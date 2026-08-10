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
 /**
 * @brief Holds development options for variance reduction
 */
struct VarianceReductionParameters
{
  enum class Strategy {
    None,
    EulerFluid,
    LocalMaxwellian,
    SpatiallyAveraged,
  };

  Strategy strategy = Strategy::None;
  double reference_number_density = 0.0;
  std::vector<double> reference_bulk_velocity{0.0,0.0,0.0};
  double reference_temperature = 0.0;
  bool specified_lf_moments = false;
  bool limit_variance_reduction = true;
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
