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

class VerletTimeIntegrator {
public:

  VerletTimeIntegrator(Discretization &discretization)
  : discretization_(discretization)
  {}

  /**
   * @brief Advance the particle and low fidelity states one timestep using the Verlet algorithm
   *
   * @param[inout] low_fidelity_states Vector of \ref LowFidelityState 's to be updated
   * @param low_fidelity_operations Vector of \ref LowFidelityOperations that form right-hand-side contributions
   * @param[inout] particle_container Particle container to be updated
   * @param particle_operations \ref ParticleOperations that form the right-hand-side contributions
   * @param[inout] field_state \ref ElectrostaticFieldState to be updated
   * @param field_operations \ref ElectrostaticFieldOperations that form right-hand-side contributions
   * @param dt Time step
   */

  void advanceTimestep(
    std::vector<LowFidelityState>& low_fidelity_states,
    const std::vector<std::unique_ptr<LowFidelityOperations>>& low_fidelity_operations,
    ParticleContainer& particle_container,
    const ParticleOperations& particle_operations,
    ElectrostaticFieldState& field_state,
    ElectrostaticFieldOperations& field_operations,
    double dt
  );
private:
  Discretization & discretization_;
};

} // namespace mfpic
