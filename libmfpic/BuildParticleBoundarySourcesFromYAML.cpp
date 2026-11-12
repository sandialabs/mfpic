#include <libmfpic/BuildParticleBoundarySourcesFromYAML.hpp>
#include <libmfpic/Errors.hpp>
#include <libmfpic/MaxwellianParticleBoundarySource.hpp>
#include <libmfpic/ParseHelpers.hpp>

#include <iostream>

namespace mfpic {

namespace {

const std::string num_particles_key = "Number of Macroparticles per Species per Timestep";

int parseNumParticlesPerTimestep(const YAML::Node& boundary_source_node) {
  const YAML::Node& num_particles_node = boundary_source_node[num_particles_key];
  if (not num_particles_node) {
    errorWithUserMessage(formatParseMessage(boundary_source_node, num_particles_key + " is required for a boundary source!"));
  }
  const int num_particles = num_particles_node.as<int>();
  if (num_particles <= 0) {
    errorWithUserMessage(formatParseMessage(num_particles_node, num_particles_key + " is nonpositive!"));
  }
  return num_particles;
}

} // namespace

std::vector<std::unique_ptr<ParticleBoundarySource>> buildParticleBoundarySourcesFromYAML(
  const YAML::Node& boundary_sources_node,
  const std::unordered_map<std::string, Species>& species_map,
  std::shared_ptr<mfem::Mesh> mesh,
  int velocity_dims)
{
  std::vector<std::unique_ptr<ParticleBoundarySource>> boundary_sources;

  if (not boundary_sources_node.IsDefined()) {
    return boundary_sources;
  }

  if (not boundary_sources_node.IsSequence()) {
    std::cout << formatParseMessage(
      boundary_sources_node, "Not a sequence, so it's ignored. Continuing.") << std::endl;
    return boundary_sources;
  }

  for (const YAML::Node& boundary_source : boundary_sources_node) {
    const std::vector<std::string> species_names = parseSpeciesNames(boundary_source["Species"], species_map);
    const int boundary_attribute = parseBoundaryAttribute(boundary_source, mesh->Dimension());
    const int num_particles_per_timestep = parseNumParticlesPerTimestep(boundary_source);

    if (boundary_source["Maxwellian"]) {
      const YAML::Node& maxwellian_node = boundary_source["Maxwellian"];
      const double flux_density = parseRequiredPositiveDouble(maxwellian_node, "Flux Density");
      const double temperature = parseRequiredPositiveDouble(maxwellian_node, "Temperature");
      for (const std::string& species_name : species_names) {
        boundary_sources.push_back(std::make_unique<MaxwellianParticleBoundarySource>(
          mesh,
          boundary_attribute,
          species_map.at(species_name),
          num_particles_per_timestep,
          flux_density,
          temperature,
          velocity_dims));
      }
    } else {
      errorWithUserMessage(formatParseMessage(boundary_source, "Missing a boundary source type."));
    }
  }

  return boundary_sources;
}

} // namespace mfpic
