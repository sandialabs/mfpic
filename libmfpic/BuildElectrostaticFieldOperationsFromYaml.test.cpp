#include <libmfpic/BuildElectrostaticFieldOperationsFromYaml.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/ElectrostaticFieldOperations.hpp>
#include <libmfpic/Pinning.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <mfem/mfem.hpp>
#include <yaml-cpp/yaml.h>

using namespace mfpic;
using namespace testing;

TEST(BuildElectrostaticFieldOperationsFromYaml, YamlWithNoBoltzmannElectronsGetsNoBoltzmannElectrons) {
  constexpr int num_elems = 10;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);
  constexpr int hgrad_order = 1;
  Discretization es_discretization(&mesh, hgrad_order);
  auto pinning = std::make_unique<Pinning>();
  const std::string fields_yaml = "";
  const YAML::Node fields_node = YAML::Load(fields_yaml);

  std::unique_ptr<ElectrostaticFieldOperations> field_operations = buildElectrostaticFieldOperationsFromYaml(
    fields_node,
    es_discretization,
    std::make_unique<Pinning>()
  );

  EXPECT_THAT(typeid(*field_operations).name(), HasSubstr("ElectrostaticFieldOperations"));
}

TEST(BuildElectrostaticFieldOperationsFromYaml, YamlWithBoltzmannElectronsGetsBoltzmannElectrons) {
  constexpr int num_elems = 10;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);
  constexpr int hgrad_order = 1;
  Discretization es_discretization(&mesh, hgrad_order);
  auto pinning = std::make_unique<Pinning>();
  const std::string fields_yaml = R"(
Boltzmann Electrons:
  Reference Number Density: 1.0e20
  Temperature: 11600.0
  )";
  const YAML::Node fields_node = YAML::Load(fields_yaml);

  std::unique_ptr<ElectrostaticFieldOperations> field_operations = buildElectrostaticFieldOperationsFromYaml(
    fields_node,
    es_discretization,
    std::make_unique<Pinning>()
  );

  EXPECT_THAT(typeid(*field_operations).name(), HasSubstr("ElectrostaticFieldOperationsWithBoltzmannElectrons"));
}
