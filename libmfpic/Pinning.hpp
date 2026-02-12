#pragma once

#include <libmfpic/DirichletBoundaryConditions.hpp>

namespace mfpic {

class Pinning : public DirichletBoundaryConditions {
public:
  /**
   * @brief Construct a new Pinning object
   *
   * @param index_of_pinned_node - which dof index to pin
   */
  Pinning(const int index_of_pinned_node = 0);

  /**
   * @brief Get a list of dof indices to pin
   *
   * @return mfem::Array<int> - list of dof indices to pin
   */
  mfem::Array<int> getDirichletBoundaryDofIndices() const;

  /**
   * @brief Apply pinning to a GridFunction
   *
   * @param solution - which GridFunction to pin
   */
  void applyBoundaryConditions(mfem::GridFunction& solution) const;

private:
  int index_of_pinned_node_;
};

}