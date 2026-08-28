#pragma once

#include <libmfpic/TimeIntegrator.hpp>

namespace mfpic {

class ElectrostaticFieldOperations;
class ElectrostaticFieldState;
class LowFidelityOperations;
class LowFidelityState;
class ParticleContainer;
class ParticleOperations;
class Discretization;

class SSPERK32TimeIntegrator : public TimeIntegrator {
public:
  SSPERK32TimeIntegrator(Discretization& discretization) : discretization_(discretization) {}

  /**
   * @brief Advance the low fidelity states one timestep using the SSPERK32 algorithm
   *  doesn't update particle data.
   *  The SSPERK32 method is a three stage second order explicit runge kutta scheme with stability on the imaginary axis.
   *
   * @param[inout] generator Random number generator.
   * @param[inout] low_fidelity_states Vector of \ref LowFidelityState 's to be updated
   * @param[inout] low_fidelity_field_states Vector of \ref ElectrostaticFieldStates to be updated
   * @param low_fidelity_operations Vector of \ref LowFidelityOperations that form right-hand-side contributions
   * @param[inout] particle_container Particle container to be updated
   * @param particle_operations \ref ParticleOperations that form the right-hand-side contributions
   * @param[in] collision_operations Vector of \ref CollisionOperations that handles collisoins.
   * @param[inout] particle_field_state Particle \ref ElectrostaticFieldState to be updated
   * @param field_operations \ref ElectrostaticFieldOperations that form right-hand-side contributions
   * @param dt Time step
   */
  virtual void advanceTimestep(
    RandomNumberGenerator& generator,
    std::vector<LowFidelityState>& low_fidelity_states,
    std::vector<ElectrostaticFieldState>& low_fidelity_field_states,
    const std::vector<std::unique_ptr<LowFidelityOperations>>& low_fidelity_operations,
    ParticleContainer& particle_container,
    ParticleOperations& particle_operations,
    const std::vector<std::unique_ptr<CollisionOperations>>& collision_operations,
    ElectrostaticFieldState& field_state,
    ElectrostaticFieldOperations& field_operations,
    double dt) const override;

private:
  Discretization& discretization_;
};

}