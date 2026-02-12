#include <libmfpic/Pinning.hpp>

namespace mfpic {

Pinning::Pinning(const int index_of_pinned_node) : index_of_pinned_node_(index_of_pinned_node) {}

mfem::Array<int> Pinning::getDirichletBoundaryDofIndices() const { return mfem::Array<int>{index_of_pinned_node_}; }

void Pinning::applyBoundaryConditions(mfem::GridFunction& solution) const {
  mfem::ConstantCoefficient zero(0.);
  mfem::Array<int> pinned_dof_indices = getDirichletBoundaryDofIndices();
  solution.ProjectCoefficient(zero, pinned_dof_indices);
}

}
