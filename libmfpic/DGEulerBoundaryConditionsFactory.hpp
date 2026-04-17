#pragma once

#include <memory>
#include <unordered_map>
#include <vector>

namespace mfem {
  class Mesh;
}

namespace YAML {
  class Node;
}

namespace mfpic {

struct DGBC;
struct Species;

enum class DGEulerBCType { REFLECTING,
                           ABSORBING };

std::unordered_map<int, DGEulerBCType> buildBoundaryAttributeToBCTypeFromYAML(
  const YAML::Node& euler_fluids_bcs,
  const int mesh_dimension);

/**
 * @brief Build boundary conditions for Euler system
 *
 * @returns Vector of vectors of \ref DGBC 's, where the outer vector is indexed by species consistent with \p species_list
 */
std::vector<std::vector<std::unique_ptr<DGBC>>> buildDGEulerBoundaryConditions(
  const std::unordered_map<int, DGEulerBCType>& boundary_attribute_to_bc_type,
  const mfem::Mesh& mesh,
  const std::vector<Species>& species_list);

}
