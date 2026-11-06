#pragma once

#include <libmfpic/Species.hpp>

#include <mfem/mfem.hpp>
#include <yaml-cpp/yaml.h>

#include <string>
#include <unordered_map>
#include <vector>

namespace mfpic {

/**
 * @brief Parse the required ``Side`` key into a boundary attribute.
 *
 * @throws If ``Side`` is missing, or is neither an inline mesh side name nor a positive integer.
 *
 * @param[in] parent_node    YAML node containing the ``Side`` key.
 * @param[in] mesh_dimension Dimension of the mesh, used to interpret inline side names.
 *
 * @returns Boundary attribute of the side.
 */
int parseBoundaryAttribute(const YAML::Node& parent_node, int mesh_dimension);

/**
 * @brief Parse a species name and check that it exists.
 *
 * @throws If the species is not in @a species_map.
 *
 * @param[in] species_node YAML scalar naming the species.
 * @param[in] species_map  Map of species names to species definitions.
 *
 * @returns Name of the species.
 */
std::string parseSpeciesName(
  const YAML::Node& species_node,
  const std::unordered_map<std::string, Species>& species_map);

/**
 * @brief Parse a sequence of species names and check that each exists.
 *
 * @throws If @a species_node is not a sequence or a species is not in @a species_map.
 *
 * @param[in] species_node YAML sequence naming the species.
 * @param[in] species_map  Map of species names to species definitions.
 *
 * @returns Names of the species.
 */
std::vector<std::string> parseSpeciesNames(
  const YAML::Node& species_node,
  const std::unordered_map<std::string, Species>& species_map);

/**
 * @brief Parse a required, positive floating-point value.
 *
 * @throws If @a key is missing or its value is nonpositive.
 *
 * @param[in] parent_node YAML node containing @a key.
 * @param[in] key         Key of the value.
 *
 * @returns The value.
 */
double parseRequiredPositiveDouble(const YAML::Node& parent_node, const std::string& key);

/**
 * @brief Parse a sequence of floating-point values into a vector.
 *
 * @throws If @a vector_node is not a sequence.
 *
 * @param[in] vector_node YAML sequence of values.
 * @param[in] name        Name of the vector, used in the error message.
 *
 * @returns Vector of the values.
 */
mfem::Vector parseVector(const YAML::Node& vector_node, const std::string& name);

/**
 * @brief Parse a sequence of three floating-point values into a vector.
 *
 * @throws If @a vector_node is not a sequence of length 3.
 *
 * @param[in] vector_node YAML sequence of values.
 * @param[in] name        Name of the vector, used in the error message.
 *
 * @returns Vector of length 3.
 */
mfem::Vector parseVector3(const YAML::Node& vector_node, const std::string& name);

} // namespace mfpic
