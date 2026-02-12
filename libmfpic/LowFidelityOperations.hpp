#pragma once

#include <libmfpic/IntegratedCharge.hpp>
#include <libmfpic/ElectromagneticFieldsEvaluator.hpp>
#include <libmfpic/LowFidelityState.hpp>

namespace mfpic {

class LowFidelityOperations {
public:
  virtual LowFidelityState accelerate(
    double dt,
    const LowFidelityState& current_state,
    const ElectromagneticFieldsEvaluator& field_provider
  ) const = 0;

  virtual LowFidelityState move(
    double dt,
    const LowFidelityState& current_state
  ) const = 0;

  virtual IntegratedCharge assembleCharge(
    const LowFidelityState& current_state
  ) const = 0;

  virtual double estimateCFL(const double & dt, const double & smallest_cell_lengthscale) const = 0;

  virtual ~LowFidelityOperations() = default;
};

} // namespace mfpic
