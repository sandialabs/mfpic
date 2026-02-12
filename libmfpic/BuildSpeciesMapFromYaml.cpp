#include <libmfpic/BuildSpeciesMapFromYaml.hpp>
#include <libmfpic/Errors.hpp>

#include <yaml-cpp/yaml.h>

#include <utility>

namespace mfpic {

std::unordered_map<std::string, Species> buildSpeciesMapFromYaml(const YAML::Node& species_nodes) {
  std::unordered_map<std::string, Species> species_map;
  if (species_nodes.IsMap()) {
    for (YAML::const_iterator it = species_nodes.begin(); it != species_nodes.end(); it++) {
      const std::string species_name = it->first.as<std::string>();
      const YAML::Node species_node = it->second;
      const double charge = species_node["Charge"].as<double>();
      const double mass = species_node["Mass"].as<double>();

      double charge_over_mass = charge / mass;
      const YAML::Node& charge_over_mass_node = species_node["Charge Over Mass"];
      if (charge_over_mass_node) {
        charge_over_mass = charge_over_mass_node.as<double>();
      }

      double specific_heat_ratio = 5. / 3.;
      const YAML::Node& specific_heat_ratio_node = species_node["Specific Heat Ratio"];
      if (specific_heat_ratio_node) {
        specific_heat_ratio = specific_heat_ratio_node.as<double>();
      }

      Species species{
        .charge = charge, .mass = mass, .charge_over_mass = charge_over_mass, .specific_heat_ratio = specific_heat_ratio};
      species_map.emplace(
        std::piecewise_construct,
        std::forward_as_tuple(species_name),
        std::forward_as_tuple(species)
      );
    }
  } else {
    errorWithUserMessage(formatParseMessage(species_nodes, "Species must be a map/object!"));
  }
  return species_map;
}

} // namespace
