#include <libmfpic/Pinning.hpp>

#include <mfem/mfem.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

TEST(Pinning, GetDirichletBoundaryDofIndicesReturnsPinnedNode) {
  constexpr int index_of_pinned_node = 7;
  Pinning pinning(index_of_pinned_node);

  mfem::Array<int> dirichlet_boundary_dof_indices = pinning.getDirichletBoundaryDofIndices();

  EXPECT_EQ(dirichlet_boundary_dof_indices.Size(), 1);
  EXPECT_EQ(dirichlet_boundary_dof_indices[0], index_of_pinned_node);
}

TEST(Pinning, ApplyBoundaryConditionsSetsPinnedNodeToZero) {
  constexpr int num_elems = 15;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);

  constexpr int hgrad_order = 1;
  mfem::H1_FECollection hgrad_finite_element_collection(hgrad_order, mesh.Dimension());
  mfem::FiniteElementSpace hgrad_finite_element_space(&mesh, &hgrad_finite_element_collection);

  mfem::GridFunction solution(&hgrad_finite_element_space);

  constexpr double constant = 1.3;
  mfem::ConstantCoefficient constant_coefficient(constant);
  solution.ProjectCoefficient(constant_coefficient);

  constexpr int index_of_pinned_node = 12;
  Pinning pinning(index_of_pinned_node);

  pinning.applyBoundaryConditions(solution);

  mfem::Vector solution_at_nodes;
  solution.GetNodalValues(solution_at_nodes);

  for (int i_node = 0; i_node < solution_at_nodes.Size(); ++i_node) {
    if (i_node == index_of_pinned_node) {
      EXPECT_EQ(solution_at_nodes[i_node], 0.);
    } else {
      EXPECT_EQ(solution_at_nodes[i_node], constant);
    }
  }
}

}