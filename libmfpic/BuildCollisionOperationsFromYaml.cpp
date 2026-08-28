#include <libmfpic/BGKRelaxationOperations.hpp>
#include <libmfpic/BuildCollisionOperationsFromYaml.hpp>
#include <libmfpic/Errors.hpp>

namespace mfpic {

std::vector<std::unique_ptr<CollisionOperations>> buildCollisionOperationsFromYaml(
  const YAML::Node& collisions_nodes,
  const std::unordered_map<std::string, Species>& species_map
) {
  std::vector<std::unique_ptr<CollisionOperations>> collision_operations;
  if (collisions_nodes.IsSequence()) {
    for (const YAML::Node& collisions_node : collisions_nodes) {
      if (collisions_node["BGK Relaxation"]) {
        const YAML::Node bgk_relaxation_node = collisions_node["BGK Relaxation"];

        const YAML::Node collision_frequency_node = bgk_relaxation_node["Collision Frequency"];
        const double collision_frequency = collision_frequency_node.as<double>();
        if (collision_frequency <= 0.0) {
          errorWithUserMessage(formatParseMessage(collision_frequency_node, "Collision Frequency must be positive!"));
        }

        const YAML::Node species_node = bgk_relaxation_node["Species"];
        const std::string species_name = species_node.as<std::string>();
        if (not species_map.contains(species_name)) {
          errorWithUserMessage(formatParseMessage(species_node, "Species was not created in the Species block!"));
        }

        collision_operations.push_back(std::make_unique<BGKRelaxationOperations>(
          collision_frequency,
          species_map.at(species_name)
        ));
      }
      else {
        errorWithUserMessage(formatParseMessage(collisions_node, "Unrecognized collision type!"));
      }
    }
  }

  return collision_operations;
}

} // namespace mfpic
