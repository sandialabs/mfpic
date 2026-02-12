#include <libmfpic/DirichletBoundaryConditionsConstant.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/MeshUtilities.hpp>

#include <mfem/fem/gridfunc.hpp>
#include <mfem/mfem.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

int getIndexOfNodeThatContainsDof(const int dof_index, const mfem::Mesh& mesh, const Discretization& discretization) {
  int index_of_node_that_contains_dof = -1;
  const int num_vertices = mesh.GetNV();
  for (int i_vertex = 0; i_vertex < num_vertices; ++i_vertex) {
    mfem::Array<int> indices_of_dofs_on_vertex;
    discretization.getFeSpace().GetVertexDofs(i_vertex, indices_of_dofs_on_vertex);
    if (indices_of_dofs_on_vertex[0] == dof_index) {
      index_of_node_that_contains_dof = i_vertex;
      break;
    }
  }
  return index_of_node_that_contains_dof;
}

double getDofLocation(const int dof_index, const mfem::Mesh& mesh, const Discretization& discretization) {
  const int index_of_node_that_contains_dof = getIndexOfNodeThatContainsDof(dof_index, mesh, discretization);
  const double* vertex_position = mesh.GetVertex(index_of_node_that_contains_dof);
  return vertex_position[0];
}

TEST(DirichletBoundaryConditionsConstant, SettingLeftBoundaryAsDirichletSelectsCorrectDofIndex) {
  constexpr int num_elems = 10;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);

  constexpr int hgrad_order = 1;
  Discretization discretization(&mesh, hgrad_order);

  std::unordered_map<std::string, int> side_name_to_boundary_attribute = getSideNameToBoundaryAttributeForInlineMeshes(
    mesh.Dimension());
  const int left_boundary_attribute = side_name_to_boundary_attribute["left"];
  std::unordered_map<int, double> boundary_attribute_to_dirichlet_value{{left_boundary_attribute, 0.}};
  DirichletBoundaryConditionsConstant constant_dirichlet_bcs(boundary_attribute_to_dirichlet_value, discretization);

  mfem::Array<int> dirichlet_boundary_dof_indices = constant_dirichlet_bcs.getDirichletBoundaryDofIndices();
  EXPECT_EQ(dirichlet_boundary_dof_indices.Size(), 1);

  const int i_dirichlet_dof = dirichlet_boundary_dof_indices[0];

  const double dof_location = getDofLocation(i_dirichlet_dof, mesh, discretization);
  EXPECT_EQ(dof_location, 0.);
}

TEST(DirichletBoundaryConditionsConstant, SettingRightBoundaryAsDirichletSelectsCorrectDofIndex) {
  constexpr int num_elems = 10;
  constexpr double length = 1.2;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems, length);

  constexpr int hgrad_order = 1;
  Discretization discretization(&mesh, hgrad_order);

  std::unordered_map<std::string, int> side_name_to_boundary_attribute = getSideNameToBoundaryAttributeForInlineMeshes(
    mesh.Dimension());
  const int right_boundary_attribute = side_name_to_boundary_attribute["right"];

  std::unordered_map<int, double> boundary_attribute_to_dirichlet_value{{right_boundary_attribute, 0.}};
  DirichletBoundaryConditionsConstant constant_dirichlet_bcs(boundary_attribute_to_dirichlet_value, discretization);

  mfem::Array<int> dirichlet_boundary_dof_indices = constant_dirichlet_bcs.getDirichletBoundaryDofIndices();
  EXPECT_EQ(dirichlet_boundary_dof_indices.Size(), 1);

  const int i_dirichlet_dof = dirichlet_boundary_dof_indices[0];
  const double dof_location = getDofLocation(i_dirichlet_dof, mesh, discretization);
  EXPECT_EQ(dof_location, length);
}

TEST(DirichletBoundaryConditionsConstant, SettingBothBoundariesAsDirichletSelectsCorrectDofIndices) {
  constexpr int num_elems = 10;
  constexpr double length = 1.5;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems, length);

  constexpr int hgrad_order = 1;
  Discretization discretization(&mesh, hgrad_order);

  std::unordered_map<std::string, int> side_name_to_boundary_attribute = getSideNameToBoundaryAttributeForInlineMeshes(
    mesh.Dimension());
  const int left_boundary_attribute = side_name_to_boundary_attribute["left"];
  const int right_boundary_attribute = side_name_to_boundary_attribute["right"];

  std::unordered_map<int, double> boundary_attribute_to_dirichlet_value{
    {left_boundary_attribute, 0.},
    {right_boundary_attribute, 0.}};
  DirichletBoundaryConditionsConstant constant_dirichlet_bcs(boundary_attribute_to_dirichlet_value, discretization);

  mfem::Array<int> dirichlet_boundary_dof_indices = constant_dirichlet_bcs.getDirichletBoundaryDofIndices();
  EXPECT_EQ(dirichlet_boundary_dof_indices.Size(), 2);

  const int i_dirichlet_dof = dirichlet_boundary_dof_indices[0];
  const int j_dirichlet_dof = dirichlet_boundary_dof_indices[1];
  const double i_dof_location = getDofLocation(i_dirichlet_dof, mesh, discretization);
  const double j_dof_location = getDofLocation(j_dirichlet_dof, mesh, discretization);
  if (i_dof_location == 0.) {
    EXPECT_EQ(j_dof_location, length);
  } else {
    EXPECT_EQ(j_dof_location, 0.);
    EXPECT_EQ(i_dof_location, length);
  }
}

TEST(DirichletBoundaryConditionsConstant, ApplyingBoundaryConditionGivesCorrectValueAtNode) {
  constexpr int num_elems = 10;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);

  constexpr int hgrad_order = 1;
  Discretization discretization(&mesh, hgrad_order);

  constexpr double left_value = 4.3;
  std::unordered_map<std::string, int> side_name_to_boundary_attribute = getSideNameToBoundaryAttributeForInlineMeshes(
    mesh.Dimension());
  const int left_boundary_attribute = side_name_to_boundary_attribute["left"];

  std::unordered_map<int, double> boundary_attribute_to_dirichlet_value{{left_boundary_attribute, left_value}};
  DirichletBoundaryConditionsConstant constant_dirichlet_bcs(boundary_attribute_to_dirichlet_value, discretization);

  mfem::GridFunction potential(&discretization.getFeSpace());
  constexpr double initial_value = 5.6;
  potential = initial_value;

  constant_dirichlet_bcs.applyBoundaryConditions(potential);

  mfem::Array<int> dirichlet_dof_indices = constant_dirichlet_bcs.getDirichletBoundaryDofIndices();
  const int i_dof_dirichlet = dirichlet_dof_indices[0];

  for (int i_dof = 0; i_dof < potential.Size(); ++i_dof) {
    if (i_dof == i_dof_dirichlet) {
      EXPECT_NE(initial_value, potential[i_dof]);
    } else {
      EXPECT_EQ(initial_value, potential[i_dof]);
    }
  }

  mfem::Vector potential_nodal_values;
  potential.GetNodalValues(potential_nodal_values);

  const int i_node = getIndexOfNodeThatContainsDof(i_dof_dirichlet, mesh, discretization);
  EXPECT_EQ(left_value, potential_nodal_values[i_node]);
}

}