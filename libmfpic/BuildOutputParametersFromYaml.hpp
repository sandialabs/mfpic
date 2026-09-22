#pragma once

#include <string>

namespace YAML {
class Node;
}

namespace mfpic {

/**
 * @brief Holds options for output dumps.
 */
struct OutputParameters {
  int output_stride = 10;

  bool output_particle_moments = true;
  bool output_particles = true;
  bool output_mesh_data = true;
  bool output_text_data = true;

  std::string particle_dump_filename = "particles.h5part";
  std::string mesh_output_folder_name = "MeshOutput";
};

/**
 * @brief Construct OutputParameters from YAML.
 *
 * @param output YAML node.
 * @return OutputParameters
 */
OutputParameters buildOutputParametersFromYAML(const YAML::Node& output);

} // namespace mfpic
