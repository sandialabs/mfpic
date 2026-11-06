#include <libmfpic/BGKRelaxationOperations.hpp>
#include <libmfpic/BuildCollisionOperationsFromYaml.hpp>
#include <libmfpic/Errors.hpp>
#include <libmfpic/ParseHelpers.hpp>

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

        const double collision_frequency = parseRequiredPositiveDouble(bgk_relaxation_node, "Collision Frequency");
        const std::string species_name = parseSpeciesName(bgk_relaxation_node["Species"], species_map);

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
