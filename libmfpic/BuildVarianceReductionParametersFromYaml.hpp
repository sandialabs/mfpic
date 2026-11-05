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
  // How elements fall back to standard PIC when the control variate is too far from the particles
  enum class CellLimiting {
    None, // Never fall back
    EveryParticle, // Fall back if any particle in the element has |1 - f_LF / f| >= 1
    Variance, // Fall back unless sum w^2 (1 - f_LF / f)^2 < sum w^2 over the element's particles
  };

  enum class Strategy {
    None,
    EulerFluid,
    UniformMaxwellian,
  };

  Strategy strategy = Strategy::None;
  double reference_number_density = 0.0;
  std::vector<double> reference_bulk_velocity{0.0,0.0,0.0};
  double reference_temperature = 0.0;
  std::vector<std::string> reference_species;
  CellLimiting cell_limiting = CellLimiting::None;
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
