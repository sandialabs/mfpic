#pragma once

#include <libmfpic/DGEulerAssembly.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/ElectromagneticFieldsEvaluator.hpp>
#include <libmfpic/IntegratedCharge.hpp>
#include <libmfpic/LowFidelityOperations.hpp>
#include <libmfpic/LowFidelityState.hpp>

#include <mfem/mfem.hpp>
#include <unordered_map>

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
   * @brief Update the state by forcing only with the electromagnetic coupling terms
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
   * @brief Update the state with fluxes and field source terms simultaneously
   *
   * @param dt Timestep
   * @param state State including dofs and species list
   * @param field_evaluator electromagnetic field evaluator
   */

  virtual LowFidelityState moveAccelerate(
    double dt,
    const LowFidelityState& state,
    const ElectromagneticFieldsEvaluator& field_evaluator
  ) const override;

  /**
   * @brief compute the rhs that would update a state with both fluxes and field source terms simultaneously
   *
   * @param state - State to be fluxed and accelerated
   * @param field_evaluator - fields that will accelerate state
   * @return LowFidelityState - rhs that can be used to update the state
   */
  virtual LowFidelityState computeRHS(
    const LowFidelityState& state,
    const ElectromagneticFieldsEvaluator& field_evaluator
  ) const override;

  /**
   * @brief Update the state due to volumetric sources
   *
   * @param dt Timestep
   * @param state State including dofs and species list
   */

  virtual LowFidelityState addVolumetricSource(
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
  * @brief Assemble charges from the fluids into the charge density per species
  *
  * @param state State including dofs and species list
  * @return IntegratedCharge - integrated charge state for a given species
  */
  virtual IntegratedCharge assembleChargePerSpecies(
    const LowFidelityState& state,
    const int ispecies
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
  * @brief Return the fluid total kinetic energy over the whole domain
  *
  * @param state Fluid state
  *
  * @return Total kinetic energy
  */
  virtual double computeTotalKineticEnergy(const LowFidelityState& state) const override;

  virtual double computeTotalCharge(const LowFidelityState& state) const override;

  /**
  * @brief Evaluates the particle distribution function for the low fidelity state at a given position and velocity.
  *
  * @param state State including dofs and species list
  * @param position Location in physical space to evaluate PDF
  * @param velocity Location in velocity space to evaluate PDF
  * @param element  Element containing position
  * @param i_species  Evaluate PDF for the given species index
  * @return PDF value
  */
  virtual double evaluateParticleDistributionFunction(
    const LowFidelityState & state,
    const mfem::Vector position,
    const mfem::Vector velocity,
    const int element,
    const int i_species) const override;

  /**
   * @brief Compute a species- and cell-wise integral required for particle number density variance reduction.
   *
   * Computes
   * \f[
   *   \int n \, \mathrm{d} \vec{x}
   * \f]
   * over each mesh cell for each species managed by @a current_state,
   * where \f$n = n(\vec{x})\f$ is the number density of the low fidelity model.
   *
   * @param[in] finite_element_space The finite element space discretizing the low fidelity model.
   * @param[in] current_state        Low fidelity state for which to compute integral.
   *
   * @returns Vector containing the integral for each species managed by @a current_state.
   */
  virtual std::unordered_map<Species,mfem::Vector> integralForVarianceReducedNumberDensity(
    mfem::FiniteElementSpace finite_element_space,
    const LowFidelityState& current_state
  ) const override;

  /**
   * @brief Compute a species- and cell-wise integral required for particle bulk velocity variance reduction.
   *
   * Computes
   * \f[
   *   \int n \vec{u} \, \mathrm{d} \vec{x}
   * \f]
   * over each mesh cell for each species managed by @a current_state,
   * where \f$n = n(\vec{x})\f$ is the number density
   * and \f$\vec{u} = \vec{u}(\vec{x})\f$ is the bulk velocity of the low fidelity model.
   *
   * @param[in] finite_element_space The finite element space discretizing the low fidelity model.
   * @param[in] current_state        Low fidelity state for which to compute integral.
   *
   * @returns DenseMatrix containing the integral for each species managed by @a current_state.
   */
  virtual std::unordered_map<Species, mfem::DenseMatrix> integralForVarianceReducedBulkVelocity(
    mfem::FiniteElementSpace finite_element_space,
    const LowFidelityState& current_state
  ) const override;

  /**
   * @brief Compute a species- and cell-wise integral required for particle temperature variance reduction.
   *
   * Computes
   * \f[
   *   \int n \left(T + \frac{m}{d k} u^2 \right) \, \mathrm{d} \vec{x}
   * \f]
   * over each mesh cell for each species managed by @a current_state,
   * where \f$n = n(\vec{x})\f$ is the number density of the low fidelity model,
   * \f$T = T(\vec{x})\f$ is the temperature of the low fidelity model,
   * \f$\vec{u} = \vec{u}(\vec{x})\f$ is the bulk velocity of the low fidelity model,
   * \f$m\f$ is the mass of the species,
   * \f$d\f$ is the number of dimensions,
   * and \f$k\f$ is Boltzmann's constant.
   *
   * @param[in] finite_element_space The finite element space discretizing the low fidelity model.
   * @param[in] current_state        Low fidelity state for which to compute integral.
   *
   * @returns Vector containing the integral for each species managed by @a current_state.
   */
  virtual std::unordered_map<Species,mfem::Vector> integralForVarianceReducedTemperature(
    mfem::FiniteElementSpace finite_element_space,
    const LowFidelityState& current_state,
    const int velocity_dims
  ) const override;

  /**
   * @brief Compute the cell-averaged number density for each species managed by a given low fidelity state.
   *
   * @param[in] current_state Low fidelity state for which to compute number density.
   *
   * @return Vector containing the cell averages for each species managed by @a current_state.
   */
  virtual std::unordered_map<Species,mfem::Vector> getCellAveragedNumberDensity(
    const LowFidelityState& current_state
  ) const override;

  /**
   * @brief Compute the cell-averaged bulk velocity for each species managed by a given low fidelity state.
   *
   * @param[in] current_state Low fidelity state for which to compute bulk velocity.
   *
   * @return DenseMatrix containing the cell averages for each species managed by @a current_state.
   */
  virtual std::unordered_map<Species,mfem::DenseMatrix> getCellAveragedBulkVelocity(
    const LowFidelityState& current_state
  ) const override;

  /**
   * @brief Compute the cell-averaged temperature for each species managed by a given low fidelity state.
   *
   * @param[in] current_state Low fidelity state for which to compute temperature.
   *
   * @return Vector containing the cell averages for each species managed by @a current_state.
   */
  virtual std::unordered_map<Species,mfem::Vector> getCellAveragedTemperature(
    const LowFidelityState& current_state
  ) const override;

private:
  Discretization & charge_discretization_;
  std::vector<std::shared_ptr<DGEulerAssembly>> dg_assemblers_;
  mutable mfem::Vector rhs_;
  mutable mfem::Vector temp_vector_;

};

} // namespace
