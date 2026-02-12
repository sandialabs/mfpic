#pragma once

#include <memory>
#include <vector>

namespace YAML {
class Node;
}

namespace mfpic {

class ParticleBoundaryFactory;

std::pair<std::vector<std::shared_ptr<ParticleBoundaryFactory>>, std::shared_ptr<ParticleBoundaryFactory>>
buildParticleBoundariesFromYaml(const YAML::Node& particles_node, int mesh_dimension);

} // namespace mfpic
