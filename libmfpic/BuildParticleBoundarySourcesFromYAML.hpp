#pragma once

#include <libmfpic/ParticleBoundarySource.hpp>
#include <libmfpic/Species.hpp>

#include <mfem/mfem.hpp>
#include <yaml-cpp/yaml.h>

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace mfpic {

/**
 * @brief Build a list of particle boundary sources from YAML.
 *
 * @param[in] boundary_sources_node YAML sequence defining the boundary sources.
 * @param[in] species_map           Map of species names to species definitions.
 * @param[in] mesh                  Mesh through whose boundary particles are injected.
 * @param[in] velocity_dims         Number of velocity dimensions.
 *
 * @returns One particle boundary source per species per entry of @a boundary_sources_node.
 */
std::vector<std::unique_ptr<ParticleBoundarySource>> buildParticleBoundarySourcesFromYAML(
  const YAML::Node& boundary_sources_node,
  const std::unordered_map<std::string, Species>& species_map,
  std::shared_ptr<mfem::Mesh> mesh,
  int velocity_dims);

} // namespace mfpic
