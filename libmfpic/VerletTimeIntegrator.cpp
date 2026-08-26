#include <libmfpic/CollisionOperations.hpp>
#include <libmfpic/ElectrostaticFieldOperations.hpp>
#include <libmfpic/ElectromagneticFieldsEvaluator.hpp>
#include <libmfpic/ElectrostaticFieldState.hpp>
#include <libmfpic/LowFidelityOperations.hpp>
#include <libmfpic/ParticleContainer.hpp>
#include <libmfpic/ParticleOperations.hpp>
#include <libmfpic/VerletTimeIntegrator.hpp>

namespace mfpic {

void VerletTimeIntegrator::advanceTimestep(
  RandomNumberGenerator& generator,
  std::vector<LowFidelityState>& low_fidelity_states,
  std::vector<ElectrostaticFieldState>& low_fidelity_field_states,
  const std::vector<std::unique_ptr<LowFidelityOperations>>& low_fidelity_operations,
  ParticleContainer& particle_container,
  ParticleOperations& particle_operations,
  const std::vector<std::unique_ptr<CollisionOperations>>& collision_operations,
  ElectrostaticFieldState& particle_field_state,
  ElectrostaticFieldOperations& field_operations,
  double dt) const
{
  IntegratedCharge particle_charge(discretization_);

  for (int i = 0; i < std::ssize(low_fidelity_operations); i++) {
    const ElectrostaticFieldState& field_state = push_lf_with_particle_fields_ ? particle_field_state : low_fidelity_field_states[i];
    IntegratedCharge low_fidelity_charge(discretization_);
    const LowFidelityOperations& operations = *low_fidelity_operations[i];
    LowFidelityState& low_fidelity_state = low_fidelity_states[i];
    low_fidelity_state = operations.accelerate(dt/2, low_fidelity_state, field_state);
    low_fidelity_state = operations.move(dt, low_fidelity_state);
    low_fidelity_charge.addCharge(operations.assembleCharge(low_fidelity_state));
    field_operations.fieldSolve(low_fidelity_field_states[i], low_fidelity_charge);
  }

  particle_container = particle_operations.accelerate(dt/2, particle_container, particle_field_state);
  particle_container = particle_operations.move(dt, particle_container);
  particle_container.cleanOutDeadParticles();
  particle_charge.addCharge(particle_operations.assembleCharge(particle_container));

  field_operations.fieldSolve(particle_field_state, particle_charge);

  for (int i = 0; i < std::ssize(low_fidelity_operations); i++) {
    const ElectrostaticFieldState& field_state = push_lf_with_particle_fields_ ? particle_field_state : low_fidelity_field_states[i];
    const LowFidelityOperations& operations = *low_fidelity_operations[i];
    LowFidelityState& low_fidelity_state = low_fidelity_states[i];
    low_fidelity_state = operations.accelerate(dt/2, low_fidelity_state, field_state);

    low_fidelity_state = operations.addVolumetricSource(dt, low_fidelity_state);
  }

  particle_container = particle_operations.accelerate(dt/2, particle_container, particle_field_state);

  for (const std::unique_ptr<CollisionOperations>& collision_operator : collision_operations) {
    collision_operator->performCollisions(dt, generator, particle_container, particle_operations);
  }
}

} // namespace mfpic
