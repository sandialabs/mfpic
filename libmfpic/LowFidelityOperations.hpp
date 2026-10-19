#pragma once

#include <libmfpic/IntegratedCharge.hpp>
#include <libmfpic/ElectromagneticFieldsEvaluator.hpp>
#include <libmfpic/LowFidelityState.hpp>

namespace mfpic {

class LowFidelityOperations {
public:
  virtual LowFidelityState accelerate(
    double dt,
    const LowFidelityState& state,
    const ElectromagneticFieldsEvaluator& field_provider
  ) const = 0;

  virtual LowFidelityState move(
    double dt,
    const LowFidelityState& state
  ) const = 0;

  virtual LowFidelityState moveAccelerate(
    double dt,
    const LowFidelityState& state,
    const ElectromagneticFieldsEvaluator& field_provider
  ) const = 0;

  virtual LowFidelityState computeRHS(
    const LowFidelityState& state,
    const ElectromagneticFieldsEvaluator& field_provider
  ) const = 0;

  virtual LowFidelityState addVolumetricSource(
    double dt,
    const LowFidelityState& state
  ) const = 0;

  virtual IntegratedCharge assembleCharge(
    const LowFidelityState& state
  ) const = 0;

  virtual IntegratedCharge assembleChargePerSpecies(
    const LowFidelityState& state,
    const int ispecies
  ) const = 0;

  virtual double evaluateParticleDistributionFunction(
    const LowFidelityState& state,
    const mfem::Vector position,
    const mfem::Vector velocity,
    const int element,
    const int i_species) const = 0;

  virtual std::unordered_map<Species, mfem::Vector> integralForVarianceReducedNumberDensity(
    mfem::FiniteElementSpace finite_element_space,
    const LowFidelityState& current_state
  ) const = 0;

  virtual std::unordered_map<Species, mfem::DenseMatrix> integralForVarianceReducedBulkVelocity(
    mfem::FiniteElementSpace finite_element_space,
    const LowFidelityState& current_state
  ) const = 0;

  virtual std::unordered_map<Species, mfem::Vector> integralForVarianceReducedTemperature(
    mfem::FiniteElementSpace finite_element_space,
    const LowFidelityState& current_state,
    const int velocity_dims
  ) const = 0;

  /**
   * @brief Compute the cell-averaged number density for each species managed by a given low fidelity state.
   *
   * @param[in] current_state Low fidelity state for which to compute number density.
   *
   * @return Vector containing the cell averages for each species managed by @a current_state.
   */
  virtual std::unordered_map<Species, mfem::Vector> getCellAveragedNumberDensity(
    const LowFidelityState& current_state
  ) const = 0;

  /**
   * @brief Compute the cell-averaged bulk velocity for each species managed by a given low fidelity state.
   *
   * @param[in] current_state Low fidelity state for which to compute bulk velocity.
   *
   * @return DenseMatrix containing the cell averages for each species managed by @a current_state.
   */
  virtual std::unordered_map<Species, mfem::DenseMatrix> getCellAveragedBulkVelocity(
    const LowFidelityState& current_state
  ) const = 0;

  /**
   * @brief Compute the cell-averaged temperature for each species managed by a given low fidelity state.
   *
   * @param[in] current_state Low fidelity state for which to compute temperature.
   *
   * @return Vector containing the cell averages for each species managed by @a current_state.
   */
  virtual std::unordered_map<Species, mfem::Vector> getCellAveragedTemperature(
    const LowFidelityState& current_state
  ) const = 0;

  virtual double estimateCFL(const double & dt, const double & smallest_cell_lengthscale) const = 0;

  virtual double computeTotalEnergy(const LowFidelityState& state) const = 0;
  virtual double computeTotalKineticEnergy(const LowFidelityState& state) const = 0;
  virtual double computeTotalCharge(const LowFidelityState& state) const = 0;

  virtual ~LowFidelityOperations() = default;
};

} // namespace mfpic
