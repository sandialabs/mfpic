#include <libmfpic/BGKRelaxationOperations.hpp>
#include <libmfpic/BuildCollisionOperationsFromYaml.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

const std::unordered_map<std::string, Species> species_map_with_electron = {
  {"electron", Species{.charge = 1.6e-19, .mass = 9.11e-31}},
};

TEST(BuildCollisionOperationsFromYaml, EmptyYAMLYieldsNoCollisionOperations) {
  std::vector<std::unique_ptr<CollisionOperations>> collision_operations = buildCollisionOperationsFromYaml(
    YAML::Node(),
    species_map_with_electron
  );

  ASSERT_TRUE(collision_operations.empty());
}

TEST(BuildCollisionOperationsFromYaml, OneCollisionOperationsCreatedWhenYAMLSequenceContainsOneObject) {
  const std::string yaml = R"(
- BGK Relaxation:
    Species: electron
    Collision Frequency: 0.1
  )";
  const YAML::Node collisions_nodes = YAML::Load(yaml);

  std::vector<std::unique_ptr<CollisionOperations>> collision_operations = buildCollisionOperationsFromYaml(
    collisions_nodes,
    species_map_with_electron
  );

  ASSERT_EQ(std::ssize(collision_operations), 1);
}

TEST(BuildCollisionOperationsFromYaml, TwoCollisionOperationsCreatedWhenYAMLSequenceContainsTwoObjects) {
  const std::string yaml = R"(
- BGK Relaxation:
    Species: electron
    Collision Frequency: 0.1
- BGK Relaxation:
    Species: electron
    Collision Frequency: 0.4
  )";
  const YAML::Node collisions_nodes = YAML::Load(yaml);

  std::vector<std::unique_ptr<CollisionOperations>> collision_operations = buildCollisionOperationsFromYaml(
    collisions_nodes,
    species_map_with_electron
  );

  ASSERT_EQ(std::ssize(collision_operations), 2);
}

} // namespace
