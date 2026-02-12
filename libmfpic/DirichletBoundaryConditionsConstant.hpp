#pragma once

#include <libmfpic/DirichletBoundaryConditions.hpp>
#include <libmfpic/Discretization.hpp>

namespace mfpic {

class DirichletBoundaryConditionsConstant : public DirichletBoundaryConditions {
public:
  /**
   * @brief Construct a Constant Dirichlet Boundary Conditions
   * 
   * @param boundary_attribute_to_dirichlet_value - map from boundary attribute to dirichlet value being set on that boundary
   * @param discretization - which finite element space will have boundary conditions applied
   *  This must be non const in order to get boundary dofs, but discretization should not be meaningfully changed.
   */
  DirichletBoundaryConditionsConstant(
    const std::unordered_map<int, double>& boundary_attribute_to_dirichlet_value,
    Discretization& discretization);

  /**
   * @brief Get a list of indices of dofs that need dirichlet boundaries applied
   * 
   * @return mfem::Array<int> list of dof indices that need dirichlet boundaries applied
   */
  mfem::Array<int> getDirichletBoundaryDofIndices() const;

  /**
   * @brief Apply boundary conditions to a grid function
   * 
   * @param solution - grid function on which to apply boundary conditions
   */
  void applyBoundaryConditions(mfem::GridFunction& solution) const;

private:
  /// list of indices that need dirichlet boundaries applied to
  mfem::Array<int> dirichlet_dof_indices_;

  /// map from boundary attribute to dirichlet value being set on that boundary
  std::unordered_map<int, double> boundary_attribute_to_dirichlet_value_;

};

}