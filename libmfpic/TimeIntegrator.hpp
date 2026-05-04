#pragma once

#include <vector>
#include <memory>

namespace mfpic {

// Forward declarations
class ElectrostaticFieldOperations;
class ElectrostaticFieldState;
class LowFidelityOperations;
class LowFidelityState;
class ParticleContainer;
class ParticleOperations;
class Discretization;

class TimeIntegrator {
public:

  virtual ~TimeIntegrator() = default;

  /**
   * @brief Advance the particle and low fidelity states one timestep
   *
   * @param[inout] low_fidelity_states Vector of \ref LowFidelityState 's to be updated
   * @param[inout] low_fidelity_field_states Vector of \ref ElectrostaticFieldStates to be updated
   * @param low_fidelity_operations Vector of \ref LowFidelityOperations that form right-hand-side contributions
   * @param[inout] particle_container Particle container to be updated
   * @param particle_operations \ref ParticleOperations that form the right-hand-side contributions
   * @param[inout] particle_field_state Particle \ref ElectrostaticFieldState to be updated
   * @param field_operations \ref ElectrostaticFieldOperations that form right-hand-side contributions
   * @param dt Time step
   */
  virtual void advanceTimestep(
    std::vector<LowFidelityState>& low_fidelity_states,
    std::vector<ElectrostaticFieldState>& low_fidelity_field_states,
    const std::vector<std::unique_ptr<LowFidelityOperations>>& low_fidelity_operations,
    ParticleContainer& particle_container,
    const ParticleOperations& particle_operations,
    ElectrostaticFieldState& field_state,
    ElectrostaticFieldOperations& field_operations,
    double dt) const = 0;
};

}