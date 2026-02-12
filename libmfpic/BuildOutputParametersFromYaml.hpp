#pragma once

#include <string>

namespace YAML {
class Node;
}

namespace mfpic {

/**
 * @brief Holds options for output dumps
 */

struct OutputParameters { 
  int output_stride = 10;
  std::string particle_dump_filename = "particles.h5part";
  std::string mesh_output_folder_name = "MeshOutput";
};

/**
 * @brief Construct OutputParameters from yaml
 * 
 * @param output - yaml node
 * @return OutputParameters 
 */
OutputParameters buildOutputParametersFromYAML(const YAML::Node& output);

}
