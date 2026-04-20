#pragma once

#include <string>

namespace mfpic {

class ParticleContainer;
class ParticleOperations;

void dumpParticles(const ParticleContainer& particles, double simulation_time, const std::string filename = "particles.h5part");

void dumpParticleMoments(ParticleOperations& particle_operations,const ParticleContainer& particles,const std::string& file_prefix,const int step,const double time); 

} // namespace mfpic
