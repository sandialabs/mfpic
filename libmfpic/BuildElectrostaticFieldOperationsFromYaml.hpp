#pragma once

#include <memory>

namespace YAML {
class Node;
}

namespace mfpic {

class DirichletBoundaryConditions;
class Discretization;
class ElectrostaticFieldOperations;

std::unique_ptr<ElectrostaticFieldOperations> buildElectrostaticFieldOperations(
  const YAML::Node& fields_node,
  Discretization& electrostatic_discretization,
  std::unique_ptr<DirichletBoundaryConditions> dirichlet_boundary_conditions
);

} // namespace mfpic
