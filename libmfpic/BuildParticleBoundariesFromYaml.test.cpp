#include <libmfpic/BuildParticleBoundariesFromYaml.hpp>
#include <libmfpic/ReflectingParticleBoundary.hpp>

#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include <memory>

namespace {

using namespace mfpic;

constexpr int one_mesh_dimension = 1;

TEST(BuildParticleBoundariesFromYaml, InvalidBoundaryConditionThrows) {
  const std::string yaml ="Boundary Conditions: {Side: 2, Type: Invalid}";
  const YAML::Node node = YAML::Load(yaml);

  EXPECT_ANY_THROW(buildParticleBoundariesFromYaml(node, one_mesh_dimension));
}

TEST(BuildParticleBoundariesFromYaml, InvalidDefaultBoundaryConditionThrows) {
  const std::string yaml = "Default Boundary Condition: Invalid";
  const YAML::Node node = YAML::Load(yaml);

  EXPECT_ANY_THROW(buildParticleBoundariesFromYaml(node, one_mesh_dimension));
}

TEST(BuildParticleBoundariesFromYaml, NonSequenceBoundaryConditionsIsIgnored) {
  const std::string yaml = "{Boundary Conditions: Invalid, Default Boundary Condition: Reflecting}";
  const YAML::Node node = YAML::Load(yaml);

  auto [boundary_factories, default_factory] = buildParticleBoundariesFromYaml(node, one_mesh_dimension);

  EXPECT_EQ(boundary_factories.size(), 0);
}

TEST(BuildParticleBoundariesFromYaml, BoundaryAttributesSetCorrectly) {
  const std::string yaml = R"(
Boundary Conditions:
  - Type: Reflecting
    Side: 1
  - Type: Reflecting
    Side: 2
Default Boundary Condition: Reflecting
  )";
  const YAML::Node node = YAML::Load(yaml);

  auto [boundary_factories, default_factory] = buildParticleBoundariesFromYaml(node, one_mesh_dimension);

  ASSERT_EQ(boundary_factories.size(), 2);
  EXPECT_EQ(boundary_factories[0]->getBoundaryAttribute(), 1);
  EXPECT_EQ(boundary_factories[1]->getBoundaryAttribute(), 2);
}

TEST(BuildParticleBoundariesFromYaml, BuiltBoundaryFactoriesHaveExpectedTypes) {
  const std::string yaml = R"(
Boundary Conditions:
  - Type: Reflecting
    Side: 1
Default Boundary Condition: Reflecting
  )";
  const YAML::Node node = YAML::Load(yaml);

  auto [boundary_factories, default_factory] = buildParticleBoundariesFromYaml(node, one_mesh_dimension);

  ASSERT_EQ(boundary_factories.size(), 1);
  EXPECT_TRUE(std::dynamic_pointer_cast<ReflectingParticleBoundaryFactory>(boundary_factories[0]));
  EXPECT_TRUE(std::dynamic_pointer_cast<ReflectingParticleBoundaryFactory>(default_factory));
}

} // namespace
