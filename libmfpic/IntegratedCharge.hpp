#pragma once

#include <libmfpic/ParticleContainer.hpp>
#include <libmfpic/Discretization.hpp>

#include <mfem/mfem.hpp>

namespace mfpic {

class IntegratedCharge {
public:

  /**
  * @brief Construct a new IntegratedCharge object
  *
  * @param discretization - Discretization object containing the finite element space
  */
  IntegratedCharge(Discretization &discretization);

  /**
  * @brief Sum another IntegratedCharge into this one.
  *
  * @param charge_to_add - IntegratedCharge to add.
  */
  void addCharge(const IntegratedCharge& charge_to_add);

  /**
   * @brief Set the whole Integrated Charge vector
   *
   * @param integrated_charge - input integrated_charge as vector
   */
  void setIntegratedCharge(const mfem::Vector& integrated_charge);

  /**
   * @brief Get the Integrated Charge as mfem::Vector
   *
   * @return mfem::Vector - integrated charge as vector
   */
  mfem::Vector getIntegratedCharge() const;

  /**
  * @brief Set a value in the Vector associated with the IntegratedCharge state
  *
  * @param dof_id - Global DOF id of the Vector to set
  * @param value - value to insert into the Vector
  */
  void setIntegratedChargeValue(int dof_id, double value);

  /**
  * @brief Add a value in the Vector associated with the IntegratedCharge state
  *
  * @param dof_id - Global DOF id of the Vector to set
  * @param value - value to add to the Vector
  */
  void addIntegratedChargeValue(int dof_id, double value);

  /**
  * @brief Set a single value for the entire Vector associated with the IntegratedCharge state
  *
  * @param value - value to specify the Vector
  */
  void setIntegratedChargeValue(double value);

  /**
  * @brief Get a value in the Vector associated with the IntegratedCharge state
  *
  * @param dof_id - Global DOF id of the Vector to set
  * @return double - value of the Vector at the index
  */
  double getIntegratedChargeValue(int dof_id);

private:

  /// vector of size number of electrostatic local dofs, that stores integrated charge
  /// integrated_charge_[i] = \int_{\Omega} \rho \phi_i dV where \phi_i is the basis function associated with dof i
  mfem::Vector integrated_charge_;
};

} // namespace mfpic
