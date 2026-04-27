#pragma once

#include <libmfpic/DGEulerAssembly.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/ElectromagneticFieldsEvaluator.hpp>
#include <libmfpic/IntegratedCharge.hpp>
#include <libmfpic/LowFidelityOperations.hpp>
#include <libmfpic/LowFidelityState.hpp>

#include <mfem/mfem.hpp>

namespace mfpic {

class DGEulerOperations : public LowFidelityOperations {

public:

  /**
  * @brief Construct a new DGEulerOperations object
  *
  * @param charge_discretization Discretization object containing the finite element space for the charge
  * @param dg_operators Operators which compute the rhs contributions for each species
  */

  DGEulerOperations(
    Discretization &charge_discretization,
    std::vector<std::shared_ptr<DGEulerAssembly>> &dg_operators
  );

  /**
   * @brief Update the state by forcing only with the source terms
   *
   * @param dt Timestep
   * @param state State including dofs and species list
   * @param field_evaluator electromagnetic field evaluator
   */

  virtual LowFidelityState accelerate(
    double dt,
    const LowFidelityState& state,
    const ElectromagneticFieldsEvaluator& field_evaluator
  ) const override;

  /**
   * @brief Update the state due to flux terms
   *
   * @param dt Timestep
   * @param state State including dofs and species list
   */

  virtual LowFidelityState move(
    double dt,
    const LowFidelityState& state
  ) const override;

  /**
  * @brief Assemble charges from the fluids into the charge density
  *
  * @param state State including dofs and species list
  * @return IntegratedCharge - integrated charge state
  */
  virtual IntegratedCharge assembleCharge(
    const LowFidelityState& state
  ) const override;

  /**
  * @brief Return the CFL based on the maximum eigenvalue of the Euler system (fluid plus acoustic speed)
  *
  * @return CFL
  */
  virtual double estimateCFL(const double & dt, const double & smallest_cell_lengthscale) const override;

/**
* @brief Return the fluid total energy over the whole domain
*
* @param state Fluid state
*
* @return Total energy
*/
  virtual double computeTotalEnergy(const LowFidelityState& state) const override;

  /**
  * @brief Evaluates the particle distribution function for the low fidelity state at a given position and velocity.
  *
  * @param state State including dofs and species list
  * @param position Location in physical space to evaluate PDF
  * @param velocity Location in velocity space to evaluate PDF
  * @param element  Element containing position
  * @param species  Evaluate PDF for the given species
  * @return PDF value
  */
  virtual double evaluateParticleDistributionFunction(
    const LowFidelityState & state,
    const mfem::Vector position,
    const mfem::Vector velocity,
    const int element,
    const Species& species) const override;

private:
  Discretization & charge_discretization_;
  std::vector<std::shared_ptr<DGEulerAssembly>> dg_assemblers_;
  mutable mfem::Vector rhs_;
  mutable mfem::Vector temp_vector_;

};

} // namespace
