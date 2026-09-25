#pragma once

#include <libmfpic/Species.hpp>

#include <string>
#include <unordered_map>

namespace YAML {
class Node;
}

namespace mfpic {

/**
 * @brief Build the map from species name to Species from the Species section of the input deck.
 *
 * @param[in] species_node   Species section of the input deck
 * @param[in] velocity_dims  Number of velocity dimensions; sets the default Specific Heat Ratio to (d + 2) / d
 *
 * @returns Map from species name to Species
 */
std::unordered_map<std::string, Species> buildSpeciesMapFromYaml(const YAML::Node& species_node, const int velocity_dims = 3);

} // namespace mfpic
