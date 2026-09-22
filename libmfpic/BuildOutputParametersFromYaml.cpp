#include <libmfpic/BuildOutputParametersFromYaml.hpp>
#include <libmfpic/Errors.hpp>

#include <yaml-cpp/yaml.h>

namespace mfpic {

OutputParameters buildOutputParametersFromYAML(const YAML::Node& output) {
  OutputParameters parameters;

  if (output["Stride"]) {
    const int stride = output["Stride"].as<int>();
    if (stride < 1) {
      errorWithUserMessage(
        formatParseMessage(output["Stride"], "Stride must be greater than 0!"));
    }
    parameters.output_stride = stride;
  }

  const auto read_bool = [&](const char* key, bool& value) {
    if (output[key]) {
      value = output[key].as<bool>();
    }
  };

  read_bool("Particle Moments", parameters.output_particle_moments);
  read_bool("Particles", parameters.output_particles);
  read_bool("Mesh Data", parameters.output_mesh_data);
  read_bool("Text Data", parameters.output_text_data);

  if (output["Particle Dump Filename"]) {
    parameters.particle_dump_filename =
      output["Particle Dump Filename"].as<std::string>();

    if (!parameters.particle_dump_filename.ends_with(".h5part")) {
      parameters.particle_dump_filename += ".h5part";
    }
  }

  if (output["Mesh Output Folder"]) {
    parameters.mesh_output_folder_name =
      output["Mesh Output Folder"].as<std::string>();
  }

  return parameters;
}

} // namespace mfpic