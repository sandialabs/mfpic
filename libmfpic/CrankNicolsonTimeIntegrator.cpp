#include <libmfpic/CrankNicolsonTimeIntegrator.hpp>

#include <libmfpic/LowFidelityOperations.hpp>

namespace mfpic {

void CrankNicolsonTimeIntegrator::advanceTimestep(
  std::vector<LowFidelityState>& low_fidelity_states,
  std::vector<ElectrostaticFieldState>& low_fidelity_field_states,
  const std::vector<std::unique_ptr<LowFidelityOperations>>& low_fidelity_operations,
  ParticleContainer& particle_container,
  const ParticleOperations& /*particle_operations*/,
  ElectrostaticFieldState& /*field_state*/,
  ElectrostaticFieldOperations& /*field_operations*/,
  double dt) const
{
  if (particle_container.numParticles() > 0) {
    errorWithUserMessage("Crank Nicolson can only be used with no PIC particles.");
  }

  for (int i_lf_model = 0; i_lf_model < std::ssize(low_fidelity_operations); i_lf_model++) {
    const LowFidelityOperations& operations = *low_fidelity_operations[i_lf_model];
    LowFidelityState& low_fidelity_state = low_fidelity_states[i_lf_model];
    ElectrostaticFieldState& low_fidelity_field_state = low_fidelity_field_states[i_lf_model];

    low_fidelity_state = operations.move(dt, low_fidelity_state);
    std::tie(low_fidelity_state, low_fidelity_field_state) = operations.plasmaOscillate(
      dt, low_fidelity_state, low_fidelity_field_state);

  }
}

}