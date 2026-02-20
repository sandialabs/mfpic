#include <libmfpic/ElectrostaticFieldOperations.hpp>
#include <libmfpic/ElectromagneticFieldsEvaluator.hpp>
#include <libmfpic/ElectrostaticFieldState.hpp>
#include <libmfpic/LowFidelityOperations.hpp>
#include <libmfpic/ParticleContainer.hpp>
#include <libmfpic/ParticleOperations.hpp>
#include <libmfpic/VerletTimeIntegrator.hpp>

namespace mfpic {

void VerletTimeIntegrator::advanceTimestep(
  std::vector<LowFidelityState>& low_fidelity_states,
  const std::vector<std::unique_ptr<LowFidelityOperations>>& low_fidelity_operations,
  ParticleContainer& particle_container,
  const ParticleOperations& particle_operations,
  ElectrostaticFieldState& field_state,
  ElectrostaticFieldOperations& field_operations,
  double dt
) {
  IntegratedCharge total_charge(discretization_);

  for (int i = 0; i < std::ssize(low_fidelity_operations); i++) {
    const LowFidelityOperations& operations = *low_fidelity_operations[i];
    LowFidelityState& low_fidelity_state = low_fidelity_states[i];
    low_fidelity_state = operations.accelerate(dt/2, low_fidelity_state, field_state);
    low_fidelity_state = operations.move(dt, low_fidelity_state);
    total_charge.addCharge(operations.assembleCharge(low_fidelity_state));
  }

  particle_container = particle_operations.accelerate(dt/2, particle_container, field_state);
  particle_container = particle_operations.move(dt, particle_container);
  total_charge.addCharge(particle_operations.assembleCharge(particle_container));

  field_operations.fieldSolve(field_state, total_charge);

  for (int i = 0; i < std::ssize(low_fidelity_operations); i++) {
    const LowFidelityOperations& operations = *low_fidelity_operations[i];
    LowFidelityState& low_fidelity_state = low_fidelity_states[i];
    low_fidelity_state = operations.accelerate(dt/2, low_fidelity_state, field_state);
  }

  particle_container = particle_operations.accelerate(dt/2, particle_container, field_state);
}

} // namespace mfpic
