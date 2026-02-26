#pragma once

#include <memory>

namespace YAML {
class Node;
}

namespace mfpic {

class DirichletBoundaryConditions;
class Discretization;
class ElectrostaticFieldOperations;

/**
  @brief Construct an electrostatic field operations object using a specification from a YAML file.

  @param[in] fields_node                   A YAML node containing the fields specifications.
  @param[in] electrostatic_discretization  Discretization used by the field solver.
  @param[in] dirichlet_boundary_conditions Dirichlet boundary conditions.

  @returns A smart pointer to electrostatic field operations.
 */
std::unique_ptr<ElectrostaticFieldOperations> buildElectrostaticFieldOperations(
  const YAML::Node& fields_node,
  Discretization& electrostatic_discretization,
  std::unique_ptr<DirichletBoundaryConditions> dirichlet_boundary_conditions
);

} // namespace mfpic
