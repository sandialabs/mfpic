#pragma once

#include <mfem.hpp>

namespace mfpic {

/**
 * @brief Parent struct for the ghost cell approach. Children's sole purpose is to override \ref 
 * setDOFsInGhost and provide a clone.
 */

struct GhostDOFSetter {
  GhostDOFSetter() {};
  virtual ~GhostDOFSetter() = default;

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
  virtual std::unique_ptr<GhostDOFSetter> clone() const = 0;
};

} // namespace mfpic
