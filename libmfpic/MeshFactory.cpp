#include <libmfpic/Errors.hpp>
#include <libmfpic/MeshFactory.hpp>

#include <mfem/mfem.hpp>

#include <yaml-cpp/yaml.h>

namespace mfpic {

std::unordered_map<std::string, mfem::Element::Type> mesh_type_string_to_element_type{
  {"line", mfem::Element::Type::SEGMENT},
  {"tri", mfem::Element::Type::TRIANGLE},
  {"quad", mfem::Element::Type::QUADRILATERAL},
  {"tet", mfem::Element::Type::TETRAHEDRON},
  {"hex", mfem::Element::Type::HEXAHEDRON}
};

std::unordered_map<std::string, int> periodic_dimension_str_to_int{
  {"x", 0},
  {"y", 1},
  {"z", 2},
};

MeshParameters buildMeshParametersFromYAML(const YAML::Node& mesh) {
  assert(mesh.IsMap());

  MeshParameters mesh_parameters;

  if (mesh["File Name"]) {
    mesh_parameters.file_name = mesh["File Name"].as<std::string>();
  } else {
    mesh_parameters.mesh_type = mesh["Type"].as<std::string>();

    for (const YAML::Node& length : mesh["Lengths"]) {
      mesh_parameters.lengths.push_back(length.as<double>());
    }

    for (const YAML::Node& num_elements_in_dim : mesh["Number of Elements"]) {
      mesh_parameters.num_elements.push_back(num_elements_in_dim.as<int>());
    }
  }
  for (const YAML::Node& periodic_dimension : mesh["Periodic Dimensions"]) {
    const std::string periodic_dimension_str = periodic_dimension.as<std::string>();
    const int periodic_dimension_int = periodic_dimension_str_to_int[periodic_dimension_str];
    mesh_parameters.periodic_dims.push_back(periodic_dimension_int);
  }

  return mesh_parameters;
}

mfem::Mesh buildMesh(const MeshParameters& mesh_parameters) {
  mfem::Mesh mesh;

  if (mesh_parameters.file_name.empty()) {
    const int mesh_dimension = mesh_parameters.lengths.size();
    constexpr bool generate_edges = true;
    mfem::Element::Type element_type = mesh_type_string_to_element_type[mesh_parameters.mesh_type];
    switch (mesh_dimension) {
      case 1:
        mesh = mfem::Mesh::MakeCartesian1D(mesh_parameters.num_elements[0], mesh_parameters.lengths[0]);
        break;
      case 2:
        mesh = mfem::Mesh::MakeCartesian2D(
          mesh_parameters.num_elements[0],
          mesh_parameters.num_elements[1],
          element_type,
          generate_edges,
          mesh_parameters.lengths[0],
          mesh_parameters.lengths[1]);
        break;
      case 3:
        mesh = mfem::Mesh::MakeCartesian3D(
          mesh_parameters.num_elements[0],
          mesh_parameters.num_elements[1],
          mesh_parameters.num_elements[2],
          element_type,
          mesh_parameters.lengths[0],
          mesh_parameters.lengths[1],
          mesh_parameters.lengths[2]);
        break;
      default:
        errorWithUserMessage("Only meshes of dimension 1, 2, or 3 are allowed.");
    }
  } else {
    mesh = mfem::Mesh(mesh_parameters.file_name);
  }
  if (not mesh_parameters.periodic_dims.empty()) {
    mfem::Vector min;
    mfem::Vector max;
    mesh.GetBoundingBox(min, max);

    std::vector<mfem::Vector> translations;
    for (const int& i_dim : mesh_parameters.periodic_dims) {
      mfem::Vector translation(mesh.Dimension());
      translation = 0.;
      translation[i_dim] = max[i_dim] - min[i_dim];
      translations.push_back(translation);
    }

    std::vector<int> periodic_vertex_mapping = mesh.CreatePeriodicVertexMapping(translations);
    mesh = mfem::Mesh::MakePeriodic(mesh, periodic_vertex_mapping);
  }

  return mesh;
}

}
