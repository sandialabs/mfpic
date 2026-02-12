#include <libmfpic/BuildOutputParametersFromYaml.hpp>

#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

namespace {

using namespace mfpic;

TEST(BuildOutputParametersFromYaml, OutputParametersSetCorrectly) {
  const std::string yaml = R"(
Particle Dump Filename: myparts.h5part
Mesh Output Folder: output_dir
Stride: 10
  )";

  const YAML::Node node = YAML::Load(yaml);

  auto output_params = buildOutputParametersFromYAML(node);

  EXPECT_EQ(output_params.mesh_output_folder_name, "output_dir");
  EXPECT_EQ(output_params.output_stride, 10);
  EXPECT_EQ(output_params.particle_dump_filename, "myparts.h5part");
}

TEST(BuildOutputParametersFromYaml, OutputParametersAddPartExtension) {
  const std::string yaml = R"(
Particle Dump Filename: myparts
  )";

  const YAML::Node node = YAML::Load(yaml);

  auto output_params = buildOutputParametersFromYAML(node);

  EXPECT_EQ(output_params.mesh_output_folder_name, "MeshOutput");
  EXPECT_EQ(output_params.output_stride, 10);
  EXPECT_EQ(output_params.particle_dump_filename, "myparts.h5part");
}

TEST(BuildOutputParametersFromYaml, PositiveStride) {
  const std::string yaml = R"(
Stride: 0
  )";

  const YAML::Node node = YAML::Load(yaml);

  EXPECT_ANY_THROW(buildOutputParametersFromYAML(node));
}

}
