#include <libmfpic/AbsorbingParticleBoundary.hpp>
#include <libmfpic/BuildParticleBoundariesFromYaml.hpp>
#include <libmfpic/Errors.hpp>
#include <libmfpic/MeshUtilities.hpp>
#include <libmfpic/ReflectingParticleBoundary.hpp>

#include <yaml-cpp/yaml.h>

#include <iostream>

namespace mfpic {

namespace {

std::shared_ptr<ParticleBoundaryFactory> createParticleBoundaryFactoryFromTypeNode(const YAML::Node& type_node) {
  std::shared_ptr<ParticleBoundaryFactory> boundary_factory;
  const std::string type = type_node.as<std::string>();
  if (type == "Reflecting") {
    boundary_factory = std::make_shared<ReflectingParticleBoundaryFactory>();
  }
  else if (type == "Absorbing") {
    boundary_factory = std::make_shared<AbsorbingParticleBoundaryFactory>();
  } else {
    errorWithUserMessage(formatParseMessage(type_node, "Invalid particle boundary type!"));
  }
  return boundary_factory;
}

} // namespace

std::pair<std::vector<std::shared_ptr<ParticleBoundaryFactory>>, std::shared_ptr<ParticleBoundaryFactory>>
buildParticleBoundariesFromYaml(const YAML::Node& particles_node, int mesh_dimension) {
  std::vector<std::shared_ptr<ParticleBoundaryFactory>> boundary_factories;
  const YAML::Node& particle_bcs = particles_node["Boundary Conditions"];
  if (particle_bcs.IsSequence()) {
    std::unordered_map<std::string, int> side_name_to_boundary_attribute = getSideNameToBoundaryAttributeForInlineMeshes(
      mesh_dimension
    );
    for (const YAML::Node& particle_bc : particle_bcs) {
      const std::string side_name = particle_bc["Side"].as<std::string>();
      int boundary_attribute;
      if (side_name_to_boundary_attribute.contains(side_name)) {
        boundary_attribute = side_name_to_boundary_attribute.at(side_name);
      } else {
        boundary_attribute = particle_bc["Side"].as<int>();
      }

      std::shared_ptr<ParticleBoundaryFactory> boundary_factory = createParticleBoundaryFactoryFromTypeNode(
        particle_bc["Type"]
      );
      boundary_factory->setBoundaryAttribute(boundary_attribute);
      boundary_factories.push_back(boundary_factory);
    }
  } else {
    std::string message = formatParseMessage(
      particle_bcs,
      "Particle boundary conditions is not a sequence, so it's ignored. Continuing."
    );
    std::cout << message << std::endl;
  }

  std::shared_ptr<ParticleBoundaryFactory> default_boundary_factory = createParticleBoundaryFactoryFromTypeNode(
    particles_node["Default Boundary Condition"]
  );

  return std::make_pair(boundary_factories, default_boundary_factory);
}

} // namespace mfpic
