#include <libmfpic/DGEulerBoundaryConditions.hpp>
#include <libmfpic/DGEulerBoundaryConditionsFactory.hpp>
#include <libmfpic/DGBC.hpp>
#include <libmfpic/DGGhostBC.hpp>
#include <libmfpic/KineticFluxBC.hpp>
#include <libmfpic/MaxwellianInfluxBC.hpp>
#include <libmfpic/MeshUtilities.hpp>

#include <gtest/gtest.h>

#include <yaml-cpp/yaml.h>

namespace {

using namespace mfpic;

TEST(DGEulerBoundaryConditionsFactory, NoEulerFluidsGivesEmptyBoundaryAttributeToBCParametersMap) {
  const int mesh_dimension = 1;
  std::string main_string("");

  YAML::Node main = YAML::Load(main_string);
  YAML::Node fluids = main["Euler Fluids"];
  YAML::Node fluid_bcs = fluids["Boundary Conditions"];

  std::unordered_map<int, std::unique_ptr<DGEulerBCParameters>> boundary_attribute_to_bc_parameters = buildBoundaryAttributeToBCParametersFromYAML(
    fluid_bcs, mesh_dimension);

  EXPECT_TRUE(boundary_attribute_to_bc_parameters.empty());
}

TEST(DGEulerBoundaryConditionsFactory, NoBoundaryConditionsGivesEmptyBoundaryAttributeToBCParametersMap) {
  const int mesh_dimension = 1;
  std::string main_string("Euler Fluids:");

  YAML::Node main = YAML::Load(main_string);
  YAML::Node fluids = main["Euler Fluids"];
  YAML::Node fluid_bcs = fluids["Boundary Conditions"];

  std::unordered_map<int, std::unique_ptr<DGEulerBCParameters>> boundary_attribute_to_bc_parameters = buildBoundaryAttributeToBCParametersFromYAML(
    fluid_bcs, mesh_dimension);

  EXPECT_TRUE(boundary_attribute_to_bc_parameters.empty());
}

TEST(DGEulerBoundaryConditionsFactory, EmptyBoundaryConditionsGivesEmptyBoundaryAttributeToBCParametersMap) {
  const int mesh_dimension = 1;
  std::string fluids_string("Boundary Conditions:");

  YAML::Node fluids = YAML::Load(fluids_string);
  YAML::Node fluid_bcs = fluids["Boundary Conditions"];

  std::unordered_map<int, std::unique_ptr<DGEulerBCParameters>> boundary_attribute_to_bc_parameters = buildBoundaryAttributeToBCParametersFromYAML(
    fluid_bcs, mesh_dimension);

  EXPECT_TRUE(boundary_attribute_to_bc_parameters.empty());
}

TEST(DGEulerBoundaryConditionsFactory, SpecifyingBCsWithSideNamesGivesCorrectBoundaryAttributeToBCParametersMap) {
  constexpr int mesh_dimension = 2;
  std::string fluids_string(
    "Boundary Conditions:\n"
    "  - Side: top\n"
    "    Type: Reflecting\n"
    "  - Side: left\n"
    "    Type: Absorbing\n"
  );

  YAML::Node fluids = YAML::Load(fluids_string);
  YAML::Node fluid_bcs = fluids["Boundary Conditions"];

  std::unordered_map<int, std::unique_ptr<DGEulerBCParameters>> boundary_attribute_to_bc_parameters = buildBoundaryAttributeToBCParametersFromYAML(
    fluid_bcs, mesh_dimension);

  EXPECT_EQ(2, boundary_attribute_to_bc_parameters.size());

  std::unordered_map<std::string, int> side_name_to_boundary_attribute = getSideNameToBoundaryAttributeForInlineMeshes(
    mesh_dimension);

  const int top_boundary_attribute = side_name_to_boundary_attribute["top"];
  EXPECT_TRUE(boundary_attribute_to_bc_parameters.contains(top_boundary_attribute));
  EXPECT_NE(nullptr, dynamic_cast<ReflectingDGEulerBCParameters*>(boundary_attribute_to_bc_parameters[top_boundary_attribute].get()));

  const int left_boundary_attribute = side_name_to_boundary_attribute["left"];
  EXPECT_TRUE(boundary_attribute_to_bc_parameters.contains(left_boundary_attribute));
  EXPECT_NE(nullptr, dynamic_cast<AbsorbingDGEulerBCParameters*>(boundary_attribute_to_bc_parameters[left_boundary_attribute].get()));
}

TEST(DGEulerBoundaryConditionsFactory, SpecifyingBCsWithBoundaryAttributesGivesCorrectBoundaryAttributeToBCParametersMap) {
  constexpr int mesh_dimension = 2;
  std::string fluids_string(
    "Boundary Conditions:\n"
    "  - Side: 1\n"
    "    Type: Reflecting\n"
    "  - Side: 2\n"
    "    Type: Absorbing\n"
  );

  YAML::Node fluids = YAML::Load(fluids_string);
  YAML::Node fluid_bcs = fluids["Boundary Conditions"];

  std::unordered_map<int, std::unique_ptr<DGEulerBCParameters>> boundary_attribute_to_bc_parameters = buildBoundaryAttributeToBCParametersFromYAML(
    fluid_bcs, mesh_dimension);

  EXPECT_EQ(2, boundary_attribute_to_bc_parameters.size());

  EXPECT_TRUE(boundary_attribute_to_bc_parameters.contains(1));
  EXPECT_NE(nullptr, dynamic_cast<ReflectingDGEulerBCParameters*>(boundary_attribute_to_bc_parameters[1].get()));

  EXPECT_TRUE(boundary_attribute_to_bc_parameters.contains(2));
  EXPECT_NE(nullptr, dynamic_cast<AbsorbingDGEulerBCParameters*>(boundary_attribute_to_bc_parameters[2].get()));
}

TEST(DGEulerBoundaryConditionsFactory, EmptyBoundaryAttributeToBCParametersGivesEmptyBoundaryConditions) {
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(10);
  std::unordered_map<int, std::unique_ptr<DGEulerBCParameters>> boundary_attribute_to_bc_parameters;
  std::vector<Species> species_list;
  std::vector<std::vector<std::unique_ptr<DGBC>>> dg_euler_bcs = buildDGEulerBoundaryConditions(boundary_attribute_to_bc_parameters, mesh, species_list);
  EXPECT_TRUE(dg_euler_bcs.empty());
}

TEST(DGEulerBoundaryConditionsFactory, BoundaryAttributeToBCParametersGivesCorrectBoundaryConditions) {
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(10);
  constexpr int boundary_attribute = 1;
  std::unordered_map<int, std::unique_ptr<DGEulerBCParameters>> boundary_attribute_to_bc_parameters;
  boundary_attribute_to_bc_parameters[boundary_attribute] = std::make_unique<ReflectingDGEulerBCParameters>();
  const Species default_species;
  std::vector<Species> species_list = {default_species, default_species};
  std::vector<std::vector<std::unique_ptr<DGBC>>> dg_euler_bcs = buildDGEulerBoundaryConditions(boundary_attribute_to_bc_parameters, mesh, species_list);
  EXPECT_EQ(2, std::ssize(dg_euler_bcs));

  for (const auto & species_bcs : dg_euler_bcs) {
    EXPECT_EQ(1, std::ssize(species_bcs));
    ASSERT_NO_THROW([[maybe_unused]] DGGhostBC& reflecting_bc = dynamic_cast<DGGhostBC&>(*species_bcs[0]));
    DGGhostBC& reflecting_bc = dynamic_cast<DGGhostBC&>(*species_bcs[0]);
    ASSERT_NO_THROW([[maybe_unused]] DGEulerReflectingBC& dof_setter = dynamic_cast<DGEulerReflectingBC&>(*(reflecting_bc.ghost_dof_setter)));
    mfem::Array<int>& boundary_attribute_has_boundary_condition = species_bcs[0]->boundary_attribute_has_boundary_condition;
    for (int i = 0; i < boundary_attribute_has_boundary_condition.Size(); ++i) {
      if (i == boundary_attribute - 1) {
        EXPECT_TRUE(boundary_attribute_has_boundary_condition[i]);
      } else {
        EXPECT_FALSE(boundary_attribute_has_boundary_condition[i]);
      }
    }
  }
}

TEST(DGEulerBoundaryConditionsFactory, MaxwellianInfluxBCParsesFluxDensityAndTemperature) {
  constexpr int mesh_dimension = 1;
  std::string fluids_string(
    "Boundary Conditions:\n"
    "  - Side: left\n"
    "    Type: Maxwellian Influx\n"
    "    Number Flux Density: 1.5e20\n"
    "    Temperature: 11600.0\n"
  );

  YAML::Node fluids = YAML::Load(fluids_string);
  std::unordered_map<int, std::unique_ptr<DGEulerBCParameters>> boundary_attribute_to_bc_parameters = buildBoundaryAttributeToBCParametersFromYAML(
    fluids["Boundary Conditions"], mesh_dimension);

  const int left_boundary_attribute = getSideNameToBoundaryAttributeForInlineMeshes(mesh_dimension)["left"];
  ASSERT_TRUE(boundary_attribute_to_bc_parameters.contains(left_boundary_attribute));
  auto* influx_parameters = dynamic_cast<MaxwellianInfluxDGEulerBCParameters*>(
    boundary_attribute_to_bc_parameters[left_boundary_attribute].get());
  ASSERT_NE(nullptr, influx_parameters);
  EXPECT_DOUBLE_EQ(1.5e20, influx_parameters->number_flux_density);
  EXPECT_DOUBLE_EQ(11600.0, influx_parameters->temperature);
}

void expectMaxwellianInfluxParseThrows(const std::string& bc_keys) {
  constexpr int mesh_dimension = 1;
  std::string fluids_string(
    "Boundary Conditions:\n"
    "  - Side: left\n"
    "    Type: Maxwellian Influx\n" + bc_keys
  );
  YAML::Node fluids = YAML::Load(fluids_string);
  EXPECT_THROW(buildBoundaryAttributeToBCParametersFromYAML(fluids["Boundary Conditions"], mesh_dimension), std::logic_error);
}

TEST(DGEulerBoundaryConditionsFactory, MaxwellianInfluxBCWithoutFluxDensityThrows) {
  expectMaxwellianInfluxParseThrows("    Temperature: 11600.0\n");
}

TEST(DGEulerBoundaryConditionsFactory, MaxwellianInfluxBCWithoutTemperatureThrows) {
  expectMaxwellianInfluxParseThrows("    Number Flux Density: 1.5e20\n");
}

TEST(DGEulerBoundaryConditionsFactory, MaxwellianInfluxBCWithNonpositiveFluxDensityThrows) {
  expectMaxwellianInfluxParseThrows("    Number Flux Density: -1.5e20\n    Temperature: 11600.0\n");
  expectMaxwellianInfluxParseThrows("    Number Flux Density:     0.0\n    Temperature: 11600.0\n");
}

TEST(DGEulerBoundaryConditionsFactory, MaxwellianInfluxBCWithNonpositiveTemperatureThrows) {
  expectMaxwellianInfluxParseThrows("    Number Flux Density: 1.5e20\n    Temperature: 0.0\n");
  expectMaxwellianInfluxParseThrows("    Number Flux Density: 1.5e20\n    Temperature: -100.0\n");
}

TEST(DGEulerBoundaryConditionsFactory, InvalidBCTypeThrows) {
  constexpr int mesh_dimension = 1;
  YAML::Node fluids = YAML::Load("Boundary Conditions:\n  - Side: left\n    Type: Sticky\n");
  EXPECT_THROW(buildBoundaryAttributeToBCParametersFromYAML(fluids["Boundary Conditions"], mesh_dimension), std::logic_error);
}

TEST(DGEulerBoundaryConditionsFactory, BCParametersMakeCorrectBCTypes) {
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(10);
  constexpr int boundary_attribute = 2;
  const Species species;

  std::unique_ptr<DGBC> reflecting_bc = ReflectingDGEulerBCParameters().makeBC(boundary_attribute, mesh, species);
  EXPECT_NE(nullptr, dynamic_cast<DGGhostBC*>(reflecting_bc.get()));

  std::unique_ptr<DGBC> absorbing_bc = AbsorbingDGEulerBCParameters().makeBC(boundary_attribute, mesh, species);
  EXPECT_NE(nullptr, dynamic_cast<KineticFluxBC*>(absorbing_bc.get()));
  EXPECT_EQ(nullptr, dynamic_cast<MaxwellianInfluxBC*>(absorbing_bc.get()));

  constexpr double number_flux_density = 1.5e20;
  constexpr double temperature = 11600.0;
  std::unique_ptr<DGBC> influx_bc =
    MaxwellianInfluxDGEulerBCParameters(number_flux_density, temperature).makeBC(boundary_attribute, mesh, species);
  auto* maxwellian_influx_bc = dynamic_cast<MaxwellianInfluxBC*>(influx_bc.get());
  ASSERT_NE(nullptr, maxwellian_influx_bc);
  EXPECT_EQ(number_flux_density, maxwellian_influx_bc->number_flux_density);
  EXPECT_EQ(temperature, maxwellian_influx_bc->temperature);
  EXPECT_TRUE(influx_bc->boundary_attribute_has_boundary_condition[boundary_attribute - 1]);
}

}
