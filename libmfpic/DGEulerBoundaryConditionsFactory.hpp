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

struct DGGhostBC;

enum class DGEulerBCType { REFLECTING };

std::unordered_map<int, DGEulerBCType> buildBoundaryAttributeToBCTypeFromYAML(
  const YAML::Node& euler_fluids_bcs,
  const int mesh_dimension);

std::vector<std::unique_ptr<DGGhostBC>> buildDGEulerBoundaryConditions(
  const std::unordered_map<int, DGEulerBCType>& boundary_attribute_to_bc_type,
  const mfem::Mesh& mesh);

}
