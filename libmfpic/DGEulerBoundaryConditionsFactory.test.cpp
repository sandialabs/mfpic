#include <libmfpic/DGEulerBoundaryConditions.hpp>
#include <libmfpic/DGEulerBoundaryConditionsFactory.hpp>
#include <libmfpic/DGGhostBC.hpp>
#include <libmfpic/MeshUtilities.hpp>

#include <gtest/gtest.h>

#include <yaml-cpp/yaml.h>

namespace {

using namespace mfpic;

TEST(DGEulerBoundaryConditionsFactory, NoEulerFluidsGivesEmptyBoundaryAttributeToBCTypeMap) {
  const int mesh_dimension = 1;
  std::string main_string("");

  YAML::Node main = YAML::Load(main_string);
  YAML::Node fluids = main["Euler Fluids"];
  YAML::Node fluid_bcs = fluids["Boundary Conditions"];

  std::unordered_map<int, DGEulerBCType> boundary_attribute_to_bc_type = buildBoundaryAttributeToBCTypeFromYAML(
    fluid_bcs, mesh_dimension);

  EXPECT_TRUE(boundary_attribute_to_bc_type.empty());
}

TEST(DGEulerBoundaryConditionsFactory, NoBoundaryConditionsGivesEmptyBoundaryAttributeToBCTypeMap) {
  const int mesh_dimension = 1;
  std::string main_string("Euler Fluids:");

  YAML::Node main = YAML::Load(main_string);
  YAML::Node fluids = main["Euler Fluids"];
  YAML::Node fluid_bcs = fluids["Boundary Conditions"];

  std::unordered_map<int, DGEulerBCType> boundary_attribute_to_bc_type = buildBoundaryAttributeToBCTypeFromYAML(
    fluid_bcs, mesh_dimension);

  EXPECT_TRUE(boundary_attribute_to_bc_type.empty());
}

TEST(DGEulerBoundaryConditionsFactory, EmptyBoundaryConditionsGivesEmptyBoundaryAttributeToBCTypeMap) {
  const int mesh_dimension = 1;
  std::string fluids_string("Boundary Conditions:");

  YAML::Node fluids = YAML::Load(fluids_string);
  YAML::Node fluid_bcs = fluids["Boundary Conditions"];

  std::unordered_map<int, DGEulerBCType> boundary_attribute_to_bc_type = buildBoundaryAttributeToBCTypeFromYAML(
    fluid_bcs, mesh_dimension);

  EXPECT_TRUE(boundary_attribute_to_bc_type.empty());
}

TEST(DGEulerBoundaryConditionsFactory, SpecifyingBCsWithSideNamesGivesCorrectBoundaryAttributeToBCTypeMap) {
  constexpr int mesh_dimension = 2;
  std::string fluids_string(
    "Boundary Conditions:\n"
    "  - Side: top\n"
    "    Type: Reflecting\n"
    "  - Side: left\n"
    "    Type: Reflecting\n"
  );

  YAML::Node fluids = YAML::Load(fluids_string);
  YAML::Node fluid_bcs = fluids["Boundary Conditions"];

  std::unordered_map<int, DGEulerBCType> boundary_attribute_to_bc_type = buildBoundaryAttributeToBCTypeFromYAML(
    fluid_bcs, mesh_dimension);

  EXPECT_EQ(2, boundary_attribute_to_bc_type.size());

  std::unordered_map<std::string, int> side_name_to_boundary_attribute = getSideNameToBoundaryAttributeForInlineMeshes(
    mesh_dimension);

  const int top_boundary_attribute = side_name_to_boundary_attribute["top"];
  EXPECT_TRUE(boundary_attribute_to_bc_type.contains(top_boundary_attribute));
  EXPECT_EQ(boundary_attribute_to_bc_type[top_boundary_attribute], DGEulerBCType::REFLECTING);

  const int left_boundary_attribute = side_name_to_boundary_attribute["left"];
  EXPECT_TRUE(boundary_attribute_to_bc_type.contains(left_boundary_attribute));
  EXPECT_EQ(boundary_attribute_to_bc_type[left_boundary_attribute], DGEulerBCType::REFLECTING);
}

TEST(DGEulerBoundaryConditionsFactory, SpecifyingBCsWithBoundaryAttributesGivesCorrectBoundaryAttributeToBCTypeMap) {
  constexpr int mesh_dimension = 2;
  std::string fluids_string(
    "Boundary Conditions:\n"
    "  - Side: 1\n"
    "    Type: Reflecting\n"
    "  - Side: 2\n"
    "    Type: Reflecting\n"
  );

  YAML::Node fluids = YAML::Load(fluids_string);
  YAML::Node fluid_bcs = fluids["Boundary Conditions"];

  std::unordered_map<int, DGEulerBCType> boundary_attribute_to_bc_type = buildBoundaryAttributeToBCTypeFromYAML(
    fluid_bcs, mesh_dimension);

  EXPECT_EQ(2, boundary_attribute_to_bc_type.size());

  EXPECT_TRUE(boundary_attribute_to_bc_type.contains(1));
  EXPECT_EQ(boundary_attribute_to_bc_type[1], DGEulerBCType::REFLECTING);

  EXPECT_TRUE(boundary_attribute_to_bc_type.contains(2));
  EXPECT_EQ(boundary_attribute_to_bc_type[2], DGEulerBCType::REFLECTING);
}

TEST(DGEulerBoundaryConditionsFactory, EmptyBoundaryAttributeToBCTypeGivesEmptyBoundaryConditions) {
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(10);
  std::unordered_map<int, DGEulerBCType> boundary_attribute_to_bc_type;
  std::vector<std::unique_ptr<DGGhostBC>> dg_euler_bcs = buildDGEulerBoundaryConditions(boundary_attribute_to_bc_type, mesh);
  EXPECT_TRUE(dg_euler_bcs.empty());
}

TEST(DGEulerBoundaryConditionsFactory, BoundaryAttributeToBCTypeGivesCorrectBoundaryConditions) {
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(10);
  constexpr int boundary_attribute = 1;
  std::unordered_map<int, DGEulerBCType> boundary_attribute_to_bc_type{{1, DGEulerBCType::REFLECTING}};
  std::vector<std::unique_ptr<DGGhostBC>> dg_euler_bcs = buildDGEulerBoundaryConditions(boundary_attribute_to_bc_type, mesh);
  EXPECT_EQ(1, std::ssize(dg_euler_bcs));

  ASSERT_NO_THROW([[maybe_unused]] DGEulerReflectingBC& reflecting_bc = dynamic_cast<DGEulerReflectingBC&>(*dg_euler_bcs[0]));

  mfem::Array<int>& boundary_attribute_has_boundary_condition = dg_euler_bcs[0]->boundary_attribute_has_boundary_condition;
  for (int i = 0; i < boundary_attribute_has_boundary_condition.Size(); ++i) {
    if (i == boundary_attribute - 1) {
      EXPECT_TRUE(boundary_attribute_has_boundary_condition[i]);
    } else {
      EXPECT_FALSE(boundary_attribute_has_boundary_condition[i]);
    }
  }

}

}