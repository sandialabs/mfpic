#include <libmfpic/Errors.hpp>
#include <libmfpic/MeshUtilities.hpp>
#include <libmfpic/ParseHelpers.hpp>

namespace mfpic {

int parseBoundaryAttribute(const YAML::Node& parent_node, int mesh_dimension) {
  const YAML::Node& side_node = parent_node["Side"];
  if (not side_node) {
    errorWithUserMessage(formatParseMessage(parent_node, "Side is required!"));
  }

  const std::unordered_map<std::string, int> side_name_to_boundary_attribute =
    getSideNameToBoundaryAttributeForInlineMeshes(mesh_dimension);
  const std::string side_name = side_node.as<std::string>();
  if (side_name_to_boundary_attribute.contains(side_name)) {
    return side_name_to_boundary_attribute.at(side_name);
  }

  int boundary_attribute = -1;
  if (not YAML::convert<int>::decode(side_node, boundary_attribute) or boundary_attribute < 1) {
    errorWithUserMessage(formatParseMessage(
      side_node, "Side must be an inline mesh side name or a positive integer boundary attribute!"));
  }
  return boundary_attribute;
}

std::string parseSpeciesName(
  const YAML::Node& species_node,
  const std::unordered_map<std::string, Species>& species_map)
{
  const std::string species_name = species_node.as<std::string>();
  if (not species_map.contains(species_name)) {
    errorWithUserMessage(formatParseMessage(species_node, "Species was not created in the Species block!"));
  }
  return species_name;
}

std::vector<std::string> parseSpeciesNames(
  const YAML::Node& species_node,
  const std::unordered_map<std::string, Species>& species_map)
{
  if (not species_node.IsSequence()) {
    errorWithUserMessage(formatParseMessage(species_node, "Species is not a sequence!"));
  }

  std::vector<std::string> species_names;
  for (const YAML::Node& species : species_node) {
    species_names.push_back(parseSpeciesName(species, species_map));
  }
  return species_names;
}

double parseRequiredPositiveDouble(const YAML::Node& parent_node, const std::string& key) {
  const YAML::Node& value_node = parent_node[key];
  if (not value_node) {
    errorWithUserMessage(formatParseMessage(parent_node, key + " is required!"));
  }
  const double value = value_node.as<double>();
  if (value <= 0.0) {
    errorWithUserMessage(formatParseMessage(value_node, key + " is nonpositive!"));
  }
  return value;
}

mfem::Vector parseVector(const YAML::Node& vector_node, const std::string& name) {
  if (not vector_node.IsSequence()) {
    errorWithUserMessage(formatParseMessage(vector_node, name + " must be a sequence."));
  }
  mfem::Vector vector(vector_node.size());
  for (int i = 0; i < vector.Size(); ++i) {
    vector[i] = vector_node[i].as<double>();
  }
  return vector;
}

mfem::Vector parseVector3(const YAML::Node& vector_node, const std::string& name) {
  if (not vector_node.IsSequence() or vector_node.size() != 3) {
    errorWithUserMessage(formatParseMessage(vector_node, name + " does not have length 3!"));
  }
  return parseVector(vector_node, name);
}

} // namespace mfpic
