#include "libmfpic/CollisionOperations.hpp"
#include <libmfpic/BuildCollisionOperationsFromYaml.hpp>

namespace mfpic {

std::vector<CollisionOperations> buildCollisionOperationsFromYaml(
  const YAML::Node& ,
  const std::unordered_map<std::string, Species>&
) {
  return std::vector<CollisionOperations>();
}

} // namespace mfpic
