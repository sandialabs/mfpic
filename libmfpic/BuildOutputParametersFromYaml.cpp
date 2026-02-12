#include <libmfpic/BuildOutputParametersFromYaml.hpp>
#include <libmfpic/Errors.hpp>

#include <yaml-cpp/yaml.h>

namespace mfpic {

OutputParameters buildOutputParametersFromYAML(const YAML::Node& output) {
  OutputParameters parameters;

  if (output["Stride"]) {
    const int stride = output["Stride"].as<int>();
    if (stride < 1) {
      errorWithUserMessage(formatParseMessage(output["Stride"], "Stride must be greater than 0!"));
    }
    parameters.output_stride = stride;
  }

  if (output["Particle Dump Filename"]) {
    const std::string filename = output["Particle Dump Filename"].as<std::string>();
    parameters.particle_dump_filename = filename;
    if (not filename.ends_with(".h5part"))
      parameters.particle_dump_filename += ".h5part";
  }

  if (output["Mesh Output Folder"]) {
    const std::string folder = output["Mesh Output Folder"].as<std::string>();
    parameters.mesh_output_folder_name = folder;
  }

  return parameters;
}

} // namespace mfpic
