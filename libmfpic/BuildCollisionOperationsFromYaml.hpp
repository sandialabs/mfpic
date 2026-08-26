#pragma once

#include <libmfpic/CollisionOperations.hpp>
#include <libmfpic/Species.hpp>

#include <yaml-cpp/yaml.h>

namespace mfpic {

std::vector<std::unique_ptr<CollisionOperations>> buildCollisionOperationsFromYaml(
  const YAML::Node& collisions_nodes,
  const std::unordered_map<std::string, Species>& species_map
);

} // namespace mfpic
