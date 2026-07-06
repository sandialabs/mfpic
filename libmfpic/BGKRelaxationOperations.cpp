#include <libmfpic/BGKRelaxationOperations.hpp>

namespace mfpic {

BGKRelaxationOperations::BGKRelaxationOperations(double collision_frequency, const Species& species_to_relax)
: collision_frequency_(collision_frequency),
  species_to_relax_(species_to_relax)
{
}

void BGKRelaxationOperations::performCollisions(
  RandomNumberGenerator& generator,
  ParticleContainer& particles,
  ParticleOperations& particle_operations
) const {

}

} // namespace mfpic
