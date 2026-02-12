#pragma once

#include <memory.h>
#include <mfem.hpp>

namespace mfpic {

 /**
 * @brief Parent struct for DG boundary conditions implemented with
 * a ghost cell approach. Children's sole purpose is to override \ref 
 * setDOFsInGhost and provide a clone.
 */
struct DGGhostBC {

  DGGhostBC() = delete;
  DGGhostBC(const int boundary_attribute, const mfem::Mesh& mesh) 
  : boundary_attribute_has_boundary_condition(mesh.bdr_attributes.Max())
  {
    boundary_attribute_has_boundary_condition = false;
    boundary_attribute_has_boundary_condition[boundary_attribute - 1] = true;
  };

  virtual ~DGGhostBC() = default;

  /**
   * @brief Sets the DOFs in the ghost cell, preparing them for the
   * Riemann solver.
   *
   * @param interior_dofs The DOFs in the real, interior boundary cell
   * @param unit_normal Outward unit normal to the boundary
   * @param[out] ghost_dofs Storage for the ghost DOFs this function must fill in
   *
   * @note The DOFs are arranged (num_dofs, num_equations).
   */
  virtual void setDOFsInGhost(const mfem::DenseMatrix & interior_dofs,
                              const mfem::Vector & unit_normal,
                              mfem::DenseMatrix & ghost_dofs) const = 0;

  /// Clone method allows for derived class to deep copy itself but return as base
  virtual std::unique_ptr<DGGhostBC> clone() const = 0;

  /// Boundary attribute numbers defining the boundaries to apply the BC
  mfem::Array<int> boundary_attribute_has_boundary_condition;

};

} // namespace
