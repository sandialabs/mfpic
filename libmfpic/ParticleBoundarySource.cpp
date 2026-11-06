#include <libmfpic/ParticleBoundarySource.hpp>

namespace mfpic {

ParticleBoundarySource::~ParticleBoundarySource() = default;

ParticleContainer ParticleBoundarySource::moveForRandomFractionOfTimestep_(
  const ParticleContainer& particles_on_boundary,
  RandomNumberGenerator& generator,
  const ParticleOperations& particle_operations,
  double dt
) {
  std::uniform_real_distribution<double> uniform_distribution(0.0, dt);
  std::vector<double> time_since_crossing;
  for (int i = 0; i < particles_on_boundary.numParticles(); i++) {
    time_since_crossing.push_back(uniform_distribution(generator));
  }

  ParticleContainer moved_particles = particle_operations.move(time_since_crossing, particles_on_boundary);
  moved_particles.cleanOutDeadParticles();
  return moved_particles;
}

} // namespace mfpic
