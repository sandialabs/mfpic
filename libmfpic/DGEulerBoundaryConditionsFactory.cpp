#include <libmfpic/DGEulerBoundaryConditions.hpp>
#include <libmfpic/DGEulerBoundaryConditionsFactory.hpp>
#include <libmfpic/DGGhostBC.hpp>
#include <libmfpic/Errors.hpp>
#include <libmfpic/MeshUtilities.hpp>

#include <memory>
#include <yaml-cpp/yaml.h>

namespace mfpic {

std::unordered_map<int, DGEulerBCType> buildBoundaryAttributeToBCTypeFromYAML(
  const YAML::Node& fluid_bcs,
  const int mesh_dimension)
{
  assert(fluid_bcs.IsSequence() or fluid_bcs.IsNull() or not fluid_bcs.IsDefined());

  std::unordered_map<std::string, int> side_name_to_boundary_attribute = getSideNameToBoundaryAttributeForInlineMeshes(
    mesh_dimension);

  std::unordered_map<int, DGEulerBCType> boundary_attribute_to_bc_type;
  for (const YAML::Node& fluid_bc : fluid_bcs) {
    const std::string side_name = fluid_bc["Side"].as<std::string>();
    int boundary_attribute;
    if (side_name_to_boundary_attribute.contains(side_name)) {
      boundary_attribute = side_name_to_boundary_attribute[side_name];
    } else {
      boundary_attribute = fluid_bc["Side"].as<int>();
      if (boundary_attribute < 1) {
        errorWithUserMessage(formatParseMessage(fluid_bc["Side"], "Side must be >= 1 if specified as boundary attribute."));
      }
    }

    const std::string bc_type_string = fluid_bc["Type"].as<std::string>();
    if (bc_type_string == "Reflecting") {
      boundary_attribute_to_bc_type[boundary_attribute] = DGEulerBCType::REFLECTING;
    } else {
      errorWithUserMessage(formatParseMessage(fluid_bc["Type"], "This boundary condition type is invalid."));
    }

  }
  return boundary_attribute_to_bc_type;
}

std::vector<std::vector<std::unique_ptr<DGBC>>> buildDGEulerBoundaryConditions(
  const std::unordered_map<int, DGEulerBCType>& boundary_attribute_to_bc_type,
  const mfem::Mesh& mesh,
  const std::vector<Species>& species_list)
{
  std::vector<std::vector<std::unique_ptr<DGBC>>> dg_euler_bcs;

  for (const Species& species : species_list) {
    std::vector<std::unique_ptr<DGBC>> species_bcs;
    for (const auto& [boundary_attribute, bc_type] : boundary_attribute_to_bc_type) {
      switch (bc_type) {
        case DGEulerBCType::REFLECTING:
          auto dof_setter = std::make_unique<DGEulerReflectingBC>();
          species_bcs.push_back(std::make_unique<DGGhostBC>(boundary_attribute, mesh, species, std::move(dof_setter)));
          break;
      }
    }
    dg_euler_bcs.push_back(std::move(species_bcs));
  }

  return dg_euler_bcs;
}

}
