#pragma once

#include <memory>
#include <unordered_map>

namespace YAML {
class Node;
}

namespace mfpic {

class DirichletBoundaryConditions;
class Discretization;

/**
 * @brief Build the map to construct dirichlet boundary conditions from input yaml
 * 
 * @param fields_bcs - yaml specifying boundary conditions, should be empty or a sequence, each item in sequence should be a map
 *  with entries Side and Value
 * @param mesh_dimension - dimension of inline mesh
 * @return std::unordered_map<int, double> map from boundary attribute to dirichlet value
 *  input to buildDirichletBoundaryConditions
 */
std::unordered_map<int, double> buildBoundaryAttributeToDirichletValueFromYAML(
  const YAML::Node& fields_bcs,
  const int mesh_dimension);

/**
 * @brief Build the dirichlet boundary conditions. If no boundary conditions are specified returns a Pinning condition.
 * 
 * @param boundary_attribute_to_dirichlet_value - map specifying constant boundary conditions
 * @param discretization - finite element discretization that is having boundary conditions applied
 * @return std::unique_ptr<DirichletBoundaryConditions> - pointer to boundary conditions
 */
std::unique_ptr<DirichletBoundaryConditions> buildDirichletBoundaryConditions(
  const std::unordered_map<int, double>& boundary_attribute_to_dirichlet_value,
  Discretization& discretization);

}