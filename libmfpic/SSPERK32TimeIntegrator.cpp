#include <libmfpic/SSPERK32TimeIntegrator.hpp>

#include <libmfpic/ElectrostaticFieldOperations.hpp>
#include <libmfpic/IntegratedCharge.hpp>
#include <libmfpic/LowFidelityOperations.hpp>

namespace mfpic {

namespace {

constexpr double a10 =  5. /  6.;
constexpr double a20 = 11. / 24.;
constexpr double a21 = 11. / 24.;

constexpr double b0 = 24. / 55.;
constexpr double b1 =  1. /  5.;
constexpr double b2 =  4. / 11.;
}

void SSPERK32TimeIntegrator::advanceTimestep(
  std::vector<LowFidelityState>& low_fidelity_states,
  std::vector<ElectrostaticFieldState>& low_fidelity_field_states,
  const std::vector<std::unique_ptr<LowFidelityOperations>>& low_fidelity_operations,
  ParticleContainer& particle_container,
  const ParticleOperations& /*particle_operations*/,
  ElectrostaticFieldState& /*field_state*/,
  ElectrostaticFieldOperations& field_operations,
  double dt) const
{
  if (particle_container.numParticles() > 0) {
    errorWithUserMessage("SSPERK32 can only be used with no PIC particles.");
  }

  for (int i_lf_model = 0; i_lf_model < std::ssize(low_fidelity_operations); i_lf_model++) {
    const LowFidelityOperations& operations = *low_fidelity_operations[i_lf_model];
    LowFidelityState& low_fidelity_state = low_fidelity_states[i_lf_model];
    const ElectrostaticFieldState& low_fidelity_field_state = low_fidelity_field_states[i_lf_model];

    LowFidelityState rhs_0 = operations.computeRHS(low_fidelity_state, low_fidelity_field_state);
    LowFidelityState low_fidelity_state_1(low_fidelity_state);
    low_fidelity_state_1.addScaledState(a10 * dt, rhs_0);

    IntegratedCharge low_fidelity_charge_1(discretization_);
    low_fidelity_charge_1.addCharge(operations.assembleCharge(low_fidelity_state_1));
    ElectrostaticFieldState low_fidelity_field_state_1 = low_fidelity_field_states[i_lf_model];
    field_operations.fieldSolve(low_fidelity_field_state_1, low_fidelity_charge_1);

    LowFidelityState rhs_1 = operations.computeRHS(low_fidelity_state_1, low_fidelity_field_state_1);

    LowFidelityState low_fidelity_state_2(low_fidelity_state);
    low_fidelity_state_2.addScaledState(a20 * dt, rhs_0);
    low_fidelity_state_2.addScaledState(a21 * dt, rhs_1);

    IntegratedCharge low_fidelity_charge_2(discretization_);
    low_fidelity_charge_2.addCharge(operations.assembleCharge(low_fidelity_state_2));
    ElectrostaticFieldState low_fidelity_field_state_2 = low_fidelity_field_states[i_lf_model];
    field_operations.fieldSolve(low_fidelity_field_state_2, low_fidelity_charge_2);

    LowFidelityState rhs_2 = operations.computeRHS(low_fidelity_state_2, low_fidelity_field_state_2);

    low_fidelity_state.addScaledState(b0 * dt, rhs_0);
    low_fidelity_state.addScaledState(b1 * dt, rhs_1);
    low_fidelity_state.addScaledState(b2 * dt, rhs_2);

    IntegratedCharge low_fidelity_charge(discretization_);
    low_fidelity_charge.addCharge(operations.assembleCharge(low_fidelity_state));
    field_operations.fieldSolve(low_fidelity_field_states[i_lf_model], low_fidelity_charge);
  }
}

}