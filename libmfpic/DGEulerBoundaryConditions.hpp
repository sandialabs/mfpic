#pragma once

#include <libmfpic/DGEulerAssembly.hpp>
#include <libmfpic/DGGhostBoundaryIntegrator.hpp>
#include <libmfpic/DGGhostBC.hpp>
#include <libmfpic/Euler.hpp>
#include <memory>

namespace mfpic {

/**
 * @brief Enforces a reflecting boundary condition by setting the ghost
 * state equal to the interior state but with normal momentum reversed.
 */

struct DGEulerReflectingBC : public DGGhostBC {

  DGEulerReflectingBC() = delete;
  DGEulerReflectingBC(const int boundary_attribute, const mfem::Mesh& mesh) :
    DGGhostBC(boundary_attribute, mesh) {};

  void setDOFsInGhost(const mfem::DenseMatrix & interior_dofs,
                      const mfem::Vector & unit_normal,
                      mfem::DenseMatrix & ghost_dofs) const override
  {
    const int num_dof = interior_dofs.NumRows();

    ghost_dofs = interior_dofs;

    for (int idof = 0; idof < num_dof; ++idof) {
      const mfem::Vector momentum {interior_dofs(idof, euler::ConservativeVariables::X_MOMENTUM_DENSITY),
                                   interior_dofs(idof, euler::ConservativeVariables::Y_MOMENTUM_DENSITY),
                                   interior_dofs(idof, euler::ConservativeVariables::Z_MOMENTUM_DENSITY) };
      const double p_dot_n = momentum * unit_normal;
      ghost_dofs(idof, euler::ConservativeVariables::X_MOMENTUM_DENSITY) -= 2. * p_dot_n * unit_normal(0);
      ghost_dofs(idof, euler::ConservativeVariables::Y_MOMENTUM_DENSITY) -= 2. * p_dot_n * unit_normal(1);
      ghost_dofs(idof, euler::ConservativeVariables::Z_MOMENTUM_DENSITY) -= 2. * p_dot_n * unit_normal(2);
    }
  }

  std::unique_ptr<DGGhostBC> clone() const override {
    return std::make_unique<DGEulerReflectingBC>(*this);
  };

};

} // namespace
