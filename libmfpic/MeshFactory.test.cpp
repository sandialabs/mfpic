#include <libmfpic/MeshFactory.hpp>

#include <mfem/mfem.hpp>

#include <yaml-cpp/yaml.h>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

TEST(MeshFactory, BuildMeshParametersFromYamlWithFileName) {
  const std::string mesh_file_name = "test.g";
  const std::string main_yaml(
    "Mesh:\n"
    "  File Name: " + mesh_file_name + "\n"
  );

  YAML::Node main = YAML::Load(main_yaml);
  YAML::Node mesh = main["Mesh"];

  const MeshParameters mesh_parameters = buildMeshParametersFromYAML(mesh);

  EXPECT_EQ(mesh_file_name, mesh_parameters.file_name);
  EXPECT_TRUE(mesh_parameters.mesh_type.empty());
  EXPECT_TRUE(mesh_parameters.lengths.empty());
  EXPECT_TRUE(mesh_parameters.num_elements.empty());
}

TEST(MeshFactory, BuildInlineMeshParametersFromYaml) {
  const std::string mesh_type = "tri";
  const std::vector<double> lengths{2.1, 3.2};
  const std::vector<int> num_elements{4, 5};
  const std::string main_yaml(
    "Mesh:\n"
    "  Type: " + mesh_type + "\n"
    "  Lengths: [" + std::to_string(lengths[0]) + ", " + std::to_string(lengths[1]) + "]\n"
    "  Number of Elements: [" + std::to_string(num_elements[0]) + ", " + std::to_string(num_elements[1]) + "]\n"
  );

  YAML::Node main = YAML::Load(main_yaml);
  YAML::Node mesh = main["Mesh"];

  const MeshParameters mesh_parameters = buildMeshParametersFromYAML(mesh);

  EXPECT_TRUE(mesh_parameters.file_name.empty());
  EXPECT_EQ(mesh_type, mesh_parameters.mesh_type);
  EXPECT_EQ(lengths, mesh_parameters.lengths);
  EXPECT_EQ(num_elements, mesh_parameters.num_elements);
}

TEST(MeshFactory, BuildInlinePeriodicDimensionFromYaml) {
  const std::string main_yaml(
    "Mesh:\n"
    "  Type: line\n"
    "  Periodic Dimensions: [x, y]\n"
  );

  YAML::Node main = YAML::Load(main_yaml);
  YAML::Node mesh = main["Mesh"];

  const MeshParameters mesh_parameters = buildMeshParametersFromYAML(mesh);
  const std::vector<int> periodic_dims = mesh_parameters.periodic_dims;

  EXPECT_EQ(2, periodic_dims.size());
  EXPECT_TRUE(std::find(periodic_dims.begin(), periodic_dims.end(), 0) != periodic_dims.end());
  EXPECT_TRUE(std::find(periodic_dims.begin(), periodic_dims.end(), 1) != periodic_dims.end());
}

TEST(MeshFactory, BuildMeshFromFileName) {
  MeshParameters mesh_parameters;
  mesh_parameters.file_name = "nonexistent_mesh";

  EXPECT_DEATH(buildMesh(mesh_parameters), "MFEM abort: Mesh file not found: nonexistent_mesh");
}

void checkMesh(const MeshParameters mesh_parameters, const int expected_num_elements) {
  mfem::Mesh mesh = buildMesh(mesh_parameters);

  const int expected_dim = mesh_parameters.lengths.size();

  EXPECT_EQ(expected_dim, mesh.Dimension());
  EXPECT_EQ(expected_num_elements, mesh.GetNE());

  mfem::Vector min;
  mfem::Vector max;
  mesh.GetBoundingBox(min, max);

  for (int i_dim = 0; i_dim < mesh.Dimension(); ++i_dim) {
    EXPECT_EQ(0., min[i_dim]);
    EXPECT_EQ(mesh_parameters.lengths[i_dim], max[i_dim]);
  }
}

TEST(MeshFactory, Build1DMeshFromParameters) {
  MeshParameters mesh_parameters;
  mesh_parameters.lengths = {5.4};
  mesh_parameters.num_elements = {10};

  checkMesh(mesh_parameters, mesh_parameters.num_elements[0]);
}

TEST(MeshFactory, Build2DQuadMeshFromParameters) {
  MeshParameters mesh_parameters;
  mesh_parameters.mesh_type = "quad";
  mesh_parameters.lengths = {5.4, 3.2};
  mesh_parameters.num_elements = {5, 6};

  const int expected_num_elements = mesh_parameters.num_elements[0] * mesh_parameters.num_elements[1];

  checkMesh(mesh_parameters, expected_num_elements);
}

TEST(MeshFactory, Build2DTriMeshFromParameters) {
  MeshParameters mesh_parameters;
  mesh_parameters.mesh_type = "tri";
  mesh_parameters.lengths = {9.3, 5.8};
  mesh_parameters.num_elements = {7, 8};

  const int expected_num_elements = mesh_parameters.num_elements[0] * mesh_parameters.num_elements[1] * 2;

  checkMesh(mesh_parameters, expected_num_elements);
}

TEST(MeshFactory, Build3DHexMeshFromParameters) {
  MeshParameters mesh_parameters;
  mesh_parameters.mesh_type = "hex";
  mesh_parameters.lengths = {9.3, 5.8, 7.6};
  mesh_parameters.num_elements = {7, 8, 9};

  const int expected_num_elements =
    mesh_parameters.num_elements[0] * mesh_parameters.num_elements[1] * mesh_parameters.num_elements[2];

  checkMesh(mesh_parameters, expected_num_elements);
}

TEST(MeshFactory, Build3DTetMeshFromParameters) {
  MeshParameters mesh_parameters;
  mesh_parameters.mesh_type = "tet";
  mesh_parameters.lengths = {9.3, 5.8, 7.6};
  mesh_parameters.num_elements = {7, 8, 9};

  const int expected_num_elements =
    mesh_parameters.num_elements[0] * mesh_parameters.num_elements[1] * mesh_parameters.num_elements[2] * 6;

  checkMesh(mesh_parameters, expected_num_elements);
}

TEST(MeshFactory, BuildPeriodicMesh1D) {
  MeshParameters mesh_parameters;
  mesh_parameters.lengths = {3.3};
  mesh_parameters.num_elements = {10};
  mesh_parameters.periodic_dims = {0};

  mfem::Mesh periodic_mesh = buildMesh(mesh_parameters);

  EXPECT_EQ(periodic_mesh.GetNV(), mesh_parameters.num_elements[0]);
}

TEST(MeshFactory, BuildPeriodicMesh2D) {
  MeshParameters mesh_parameters;
  mesh_parameters.mesh_type = "quad";
  mesh_parameters.lengths = {4.5, 6.7};
  mesh_parameters.num_elements = {10, 20};
  mesh_parameters.periodic_dims = {1};

  mfem::Mesh periodic_mesh = buildMesh(mesh_parameters);

  const int expected_num_vertices = (mesh_parameters.num_elements[0] + 1) * mesh_parameters.num_elements[1];
  EXPECT_EQ(expected_num_vertices, periodic_mesh.GetNV());
}

TEST(MeshFactory, BuildPeriodicMesh3D) {
  MeshParameters mesh_parameters;
  mesh_parameters.mesh_type = "hex";
  mesh_parameters.lengths = {1.3, 4.6, 7.9};
  mesh_parameters.num_elements = {10, 20, 30};
  mesh_parameters.periodic_dims = {0, 2};

  mfem::Mesh periodic_mesh = buildMesh(mesh_parameters);

  const int expected_num_vertices =
    mesh_parameters.num_elements[0] * (mesh_parameters.num_elements[1] + 1) * mesh_parameters.num_elements[2];
  EXPECT_EQ(expected_num_vertices, periodic_mesh.GetNV());
}

}