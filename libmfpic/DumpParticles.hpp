#pragma once

#include <string>

namespace mfem {
    class Mesh;
} // namespace mfem

namespace mfpic {

class ParticleContainer;
class ParticleOperations;
class LowFidelityState;
class LowFidelityOperations;

void dumpParticles(const ParticleContainer& particles, double simulation_time, const std::string filename = "particles.h5part");

void dumpParticleMoments(ParticleOperations& particle_operations,const ParticleContainer& particles,const std::string& file_prefix,const int step,const double time);

/**
 * @brief Append variance reduced particle moments to <file_prefix>_<species>.csv
 *
 * The low-fidelity state and operations may be null for Strategy: Uniform Maxwellian.
 */
void dumpVarianceReducedParticleMoments(
    ParticleOperations& particle_operations,
    const ParticleContainer& particles,
    const LowFidelityState* low_fidelity_state,
    const LowFidelityOperations* low_fidelity_operations,
    const std::string& file_prefix,
    const int step,
    const double time);

void dumpLowFidelityMoments(
    const LowFidelityState& low_fidelity_state,
    const LowFidelityOperations& low_fidelity_operations,
    mfem::Mesh& mesh,
    const std::string& file_prefix,
    const int step,
    const double time);

} // namespace mfpic
