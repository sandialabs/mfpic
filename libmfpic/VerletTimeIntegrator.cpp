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
  auto& message = std::cout;
  message << "advanceTimeStep" << std::endl;
  IntegratedCharge total_charge(discretization_);

  message << "low fidelity accelerate and move" << std::endl;
  for (int i = 0; i < std::ssize(low_fidelity_operations); i++) {
    const LowFidelityOperations& operations = *low_fidelity_operations[i];
    LowFidelityState& low_fidelity_state = low_fidelity_states[i];
    message << "low fidelity accelerate" << std::endl;
    low_fidelity_state = operations.accelerate(dt/2, low_fidelity_state, field_state);
    message << "low fidelity move" << std::endl;
    low_fidelity_state = operations.move(dt, low_fidelity_state);
    message << "add charge" << std::endl;
    total_charge.addCharge(operations.assembleCharge(low_fidelity_state));
  }

  message << "particle accelerate" << std::endl;
  particle_container = particle_operations.accelerate(dt/2, particle_container, field_state);
  message << "particle move" << std::endl;
  particle_container = particle_operations.move(dt, particle_container);
  message << "integrated charge" << std::endl;
  total_charge.addCharge(particle_operations.assembleCharge(particle_container));

  message << "field solve" << std::endl;
  field_operations.fieldSolve(field_state, total_charge);

  message << "low fidelity accelerate" << std::endl;
  for (int i = 0; i < std::ssize(low_fidelity_operations); i++) {
    const LowFidelityOperations& operations = *low_fidelity_operations[i];
    LowFidelityState& low_fidelity_state = low_fidelity_states[i];
    low_fidelity_state = operations.accelerate(dt/2, low_fidelity_state, field_state);
  }

  message << "particle accelerate" << std::endl;
  particle_container = particle_operations.accelerate(dt/2, particle_container, field_state);
}

} // namespace mfpic
