#include <libmfpic/DirichletBoundaryConditionsFactory.hpp>
#include <libmfpic/DirichletBoundaryConditionsConstant.hpp>
#include <libmfpic/ParseHelpers.hpp>
#include <libmfpic/Pinning.hpp>

#include <yaml-cpp/yaml.h>

namespace mfpic {

std::unordered_map<int, double> buildBoundaryAttributeToDirichletValueFromYAML(
  const YAML::Node& field_bcs,
  const int mesh_dimension)
{
  assert(field_bcs.IsSequence() or field_bcs.IsNull() or not field_bcs.IsDefined());

  std::unordered_map<int, double> boundary_attribute_to_dirichlet_value;
  for (const YAML::Node& field_bc : field_bcs) {
    const int boundary_attribute = parseBoundaryAttribute(field_bc, mesh_dimension);
    const double dirichlet_value = field_bc["Value"].as<double>();
    boundary_attribute_to_dirichlet_value[boundary_attribute] = dirichlet_value;
  }

  return boundary_attribute_to_dirichlet_value;
}

std::unique_ptr<DirichletBoundaryConditions> buildDirichletBoundaryConditions(
  const std::unordered_map<int, double>& boundary_attribute_to_dirichlet_value,
  Discretization& discretization)
{
  std::unique_ptr<DirichletBoundaryConditions> dirichlet_bcs;
  if (boundary_attribute_to_dirichlet_value.empty()) {
    dirichlet_bcs = std::make_unique<Pinning>();
  } else {
    dirichlet_bcs = std::make_unique<DirichletBoundaryConditionsConstant>(boundary_attribute_to_dirichlet_value, discretization);
  }
  return dirichlet_bcs;
}

}
