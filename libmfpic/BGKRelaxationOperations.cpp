#include <libmfpic/BGKRelaxationOperations.hpp>

namespace mfpic {

BGKRelaxationOperations::BGKRelaxationOperations(double collision_frequency) : collision_frequency_(collision_frequency)
{
}

void BGKRelaxationOperations::performCollisions(
  RandomNumberGenerator& ,
  ParticleContainer& ,
  LowFidelityState&
) const {

}

} // namespace mfpic
