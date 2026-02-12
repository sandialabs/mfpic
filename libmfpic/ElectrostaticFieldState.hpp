#pragma once

#include <libmfpic/Discretization.hpp>
#include <libmfpic/ElectromagneticFieldsEvaluator.hpp>

namespace mfpic {

class ElectrostaticFieldState : public ElectromagneticFieldsEvaluator {

public:
  /**
   * @brief Construct a new Electrostatic Field State object
   *
   * @param electrostatic_discretization - the discretization for the es potential
   *  This object is non const because constructing GridFunctions requires a non const FiniteElementSpace
   *  The grid function should not meaningfully change the FiniteElementSpace
   */
  ElectrostaticFieldState(Discretization& electrostatic_discretizaton);

  /**
   * @brief Set the electrostatic potential
   * @note This only copies the basic vector data in the input potential to the stored potential.
   *  It does not update the finite element space associated with the stored potential
   *
   * @param potential - grid function of the electrostatic potential
   */
  void setPotential(const mfem::GridFunction& potential) { potential_ = potential; }

  /**
   * @brief Get the potential of the electrostatic field
   *
   * @return mfem::GridFunction - specifies potential over the mesh
   */
  template <class Self>
  auto&& getPotential(this Self&& self) { return self.potential_; }

  /**
   * @brief get the electric field that corresponds to this field state at a point in the mesh
   *
   * @param position - point in the mesh to get the electric field
   * @param element_index - index of the mesh element that contains the point
   * @return mfem::Vector - 3D vector that contains the electric field
   */
  mfem::Vector getEFieldAt(const mfem::Vector& position, const int element_index) const;

  /**
   * @brief get the B field that corresponds to this field state at a point in the mesh
   *
   * @param position - point in the mesh to get the B field
   * @param element_index - index of the mesh element that contains the point
   * @return mfem::Vector - 3D vector that contains the B field
   */
  mfem::Vector getBFieldAt(const mfem::Vector& position, const int element_index) const;

private:
  /// value of electrostatic potential on mesh
  mfem::GridFunction potential_;

};

} // namespace mfpic
