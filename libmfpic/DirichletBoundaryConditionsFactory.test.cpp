#include <libmfpic/DirichletBoundaryConditions.hpp>
#include <libmfpic/DirichletBoundaryConditionsConstant.hpp>
#include <libmfpic/DirichletBoundaryConditionsFactory.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/MeshUtilities.hpp>
#include <libmfpic/Pinning.hpp>

#include <mfem/mfem.hpp>

#include <yaml-cpp/yaml.h>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

TEST(DirichletBoundaryConditionsFactory, NoFieldsGivesEmptyBoundaryAttributeToDirichletValueMap) {
  const int mesh_dimension = 1;
  std::string main_string("");

  YAML::Node main = YAML::Load(main_string);
  YAML::Node fields = main["Fields"];
  YAML::Node field_bcs = fields["Boundary Conditions"];

  std::unordered_map<int, double> boundary_attribute_to_dirichlet_value = buildBoundaryAttributeToDirichletValueFromYAML(
    field_bcs, mesh_dimension);

  EXPECT_TRUE(boundary_attribute_to_dirichlet_value.empty());
}

TEST(DirichletBoundaryConditionsFactory, NoBoundaryConditionsGivesEmptyBoundaryAttributeToDirichletValueMap) {
  const int mesh_dimension = 1;
  std::string main_string("Fields:");

  YAML::Node main = YAML::Load(main_string);
  YAML::Node fields = main["Fields"];
  YAML::Node field_bcs = fields["Boundary Conditions"];

  std::unordered_map<int, double> boundary_attribute_to_dirichlet_value = buildBoundaryAttributeToDirichletValueFromYAML(
    field_bcs, mesh_dimension);

  EXPECT_TRUE(boundary_attribute_to_dirichlet_value.empty());
}

TEST(DirichletBoundaryConditionsFactory, EmptyBoundaryConditionsGivesEmptyBoundaryAttributeToDirichletValueMap) {
  const int mesh_dimension = 1;
  std::string fields_string("Boundary Conditions:");

  YAML::Node fields = YAML::Load(fields_string);
  YAML::Node field_bcs = fields["Boundary Conditions"];

  std::unordered_map<int, double> boundary_attribute_to_dirichlet_value = buildBoundaryAttributeToDirichletValueFromYAML(
    field_bcs, mesh_dimension);

  EXPECT_TRUE(boundary_attribute_to_dirichlet_value.empty());
}

TEST(DirichletBoundaryConditionsFactory, SpecifyingBCsWithSideNamesGiveCorrectBoundaryAttributeToDirichletValueMap) {
  constexpr int mesh_dimension = 1;
  constexpr double left_dirichlet_value = 1.2;
  constexpr double right_dirichlet_value = 2.1;
  std::string fields_string(
    "Boundary Conditions:\n"
    "  - Side: left\n"
    "    Value: " + std::to_string(left_dirichlet_value) + "\n"
    "  - Side: right\n"
    "    Value: " + std::to_string(right_dirichlet_value) + "\n"
  );

  YAML::Node fields = YAML::Load(fields_string);
  YAML::Node field_bcs = fields["Boundary Conditions"];

  std::unordered_map<int, double> boundary_attribute_to_dirichlet_value = buildBoundaryAttributeToDirichletValueFromYAML(
    field_bcs, mesh_dimension);

  EXPECT_EQ(2, boundary_attribute_to_dirichlet_value.size());

  std::unordered_map<std::string, int> side_name_to_boundary_attribute = getSideNameToBoundaryAttributeForInlineMeshes(
    mesh_dimension);

  const int left_boundary_attribute = side_name_to_boundary_attribute["left"];
  EXPECT_TRUE(boundary_attribute_to_dirichlet_value.contains(left_boundary_attribute));
  EXPECT_EQ(boundary_attribute_to_dirichlet_value[left_boundary_attribute], left_dirichlet_value);

  const int right_boundary_attribute = side_name_to_boundary_attribute["right"];
  EXPECT_TRUE(boundary_attribute_to_dirichlet_value.contains(right_boundary_attribute));
  EXPECT_EQ(boundary_attribute_to_dirichlet_value[right_boundary_attribute], right_dirichlet_value);
}

TEST(DirichletBoundaryConditionsFactory, SpecifyingBCsWithAttributesGiveCorrectBoundaryAttributeToDirichletValueMap) {
  constexpr int mesh_dimension = 1;
  constexpr double dirichlet_value_2 = 1.2;
  constexpr double dirichlet_value_3 = 2.1;
  std::string fields_string(
    "Boundary Conditions:\n"
    "  - Side: 2\n"
    "    Value: " + std::to_string(dirichlet_value_2) + "\n"
    "  - Side: 3\n"
    "    Value: " + std::to_string(dirichlet_value_3) + "\n"
  );

  YAML::Node fields = YAML::Load(fields_string);
  YAML::Node field_bcs = fields["Boundary Conditions"];

  std::unordered_map<int, double> boundary_attribute_to_dirichlet_value = buildBoundaryAttributeToDirichletValueFromYAML(
    field_bcs, mesh_dimension);

  EXPECT_EQ(2, boundary_attribute_to_dirichlet_value.size());

  EXPECT_TRUE(boundary_attribute_to_dirichlet_value.contains(2));
  EXPECT_EQ(boundary_attribute_to_dirichlet_value[2], dirichlet_value_2);

  EXPECT_TRUE(boundary_attribute_to_dirichlet_value.contains(3));
  EXPECT_EQ(boundary_attribute_to_dirichlet_value[3], dirichlet_value_3);
}

TEST(DirichletBoundaryConditionsFactory, NoConstantBoundaryConditionsReturnsPinning) {
  constexpr int num_elems = 20;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);

  constexpr int basis_order = 1;
  Discretization electrostatic_discretization(&mesh, basis_order);
  std::unordered_map<int, double> boundary_attribute_to_dirichlet_value;

  std::unique_ptr<DirichletBoundaryConditions> dirichlet_bcs = buildDirichletBoundaryConditions(
    boundary_attribute_to_dirichlet_value, electrostatic_discretization);

  ASSERT_NO_THROW([[maybe_unused]] Pinning& pinning = dynamic_cast<Pinning&>(*dirichlet_bcs));
}

TEST(DirichletBoundaryConditionsFactory, ConstantBoundaryConditionsReturnDirichletBoundaryConditionsConstant) {
  constexpr int num_elems = 20;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);

  constexpr int basis_order = 1;
  Discretization electrostatic_discretization(&mesh, basis_order);
  std::unordered_map<int, double> boundary_attribute_to_dirichlet_value;
  constexpr double dirichlet_value = 2.1;
  boundary_attribute_to_dirichlet_value[1] = dirichlet_value;

  std::unique_ptr<DirichletBoundaryConditions> dirichlet_bcs = buildDirichletBoundaryConditions(
    boundary_attribute_to_dirichlet_value, electrostatic_discretization);

  ASSERT_NO_THROW([[maybe_unused]] DirichletBoundaryConditionsConstant& constant_bcs = dynamic_cast<DirichletBoundaryConditionsConstant&>(*dirichlet_bcs));
}

}