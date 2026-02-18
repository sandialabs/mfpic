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
  particle_container.cleanOutDeadParticles();
  double ion_weight = 0.0;
  int num_ions = 0;
  int num_electons = 0;
  for (const Particle& particle : particle_container) {
    if (particle.species.charge > 0) {
      ion_weight += particle.weight;
      num_ions += 1;
    }
    else {
      num_electons += 1;
    }
  }
  std::cout << "LTM weight is " << ion_weight << std::endl;
  std::cout << "LTM num ions " << num_ions << " num electrons " << num_electons << std::endl;
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
