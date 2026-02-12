#pragma once

#include <libmfpic/Species.hpp>

#include <string>
#include <unordered_map>

namespace YAML {
class Node;
}

namespace mfpic {

std::unordered_map<std::string, Species> buildSpeciesMapFromYaml(const YAML::Node& species_node);

} // namespace mfpic
