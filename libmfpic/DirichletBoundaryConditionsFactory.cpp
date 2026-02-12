#include <libmfpic/DirichletBoundaryConditionsFactory.hpp>
#include <libmfpic/DirichletBoundaryConditionsConstant.hpp>
#include <libmfpic/MeshUtilities.hpp>
#include <libmfpic/Pinning.hpp>

#include <yaml-cpp/yaml.h>

namespace mfpic {

std::unordered_map<int, double> buildBoundaryAttributeToDirichletValueFromYAML(
  const YAML::Node& field_bcs,
  const int mesh_dimension)
{
  assert(field_bcs.IsSequence() or field_bcs.IsNull() or not field_bcs.IsDefined());

  std::unordered_map<std::string, int> side_name_to_boundary_attribute = getSideNameToBoundaryAttributeForInlineMeshes(
    mesh_dimension);

  std::unordered_map<int, double> boundary_attribute_to_dirichlet_value;
  for (const YAML::Node& field_bc : field_bcs) {
    const std::string side_name = field_bc["Side"].as<std::string>();
    int boundary_attribute;
    if (side_name_to_boundary_attribute.contains(side_name)) {
      boundary_attribute = side_name_to_boundary_attribute[side_name];
    } else {
      boundary_attribute = field_bc["Side"].as<int>();
    }

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