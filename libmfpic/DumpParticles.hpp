#pragma once

#include <string>

namespace mfpic {

class ParticleContainer;

void dumpParticles(const ParticleContainer& particles, double simulation_time, const std::string filename = "particles.h5part");

} // namespace mfpic
