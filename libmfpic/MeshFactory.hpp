#pragma once

#include <string>
#include <vector>

namespace mfem {
class Mesh;
}

namespace YAML {
class Node;
}

namespace mfpic {

struct MeshParameters {
  std::string file_name{};

  std::string mesh_type;
  std::vector<double> lengths;
  std::vector<int> num_elements;

  std::vector<int> periodic_dims{};
};

/**
 * @brief Build mesh parameters from user input as YAML
 * 
 * @param mesh - YAML specifying mesh, should be a Map that either specifies mesh with file through key "File Name" or specifies
 *  an inline mesh with keys "Type", "Lengths", and "Number of Elements"
 *  Periodicity can be specified through the key "Periodic Dimensions"
 * @return MeshParameters - parameters to create a Mesh
 */
MeshParameters buildMeshParametersFromYAML(const YAML::Node& mesh);

/**
 * @brief build a mesh from parameters
 * 
 * @param mesh_parameters - parameters defining the mesh
 * @return mfem::Mesh - mesh
 */
mfem::Mesh buildMesh(const MeshParameters& mesh_parameters);

}