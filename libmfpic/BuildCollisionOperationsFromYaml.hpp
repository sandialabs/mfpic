#pragma once

#include <libmfpic/CollisionOperations.hpp>
#include <libmfpic/Species.hpp>

#include <yaml-cpp/yaml.h>

namespace mfpic {

/**
 * @brief Build a set of collision operations objects from a YAML node.
 *
 * @param[in] collisions_node YAML node containing a sequence with collision parameters.
 * @param[in] species_map     A map of species name to species definitions.
 *
 * @return An array of collision operations.
 */
std::vector<std::unique_ptr<CollisionOperations>> buildCollisionOperationsFromYaml(
  const YAML::Node& collisions_nodes,
  const std::unordered_map<std::string, Species>& species_map
);

} // namespace mfpic
