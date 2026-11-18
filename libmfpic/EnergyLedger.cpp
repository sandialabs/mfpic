#include <libmfpic/EnergyLedger.hpp>

#include <algorithm>
#include <iomanip>
#include <iostream>

namespace mfpic {

EnergyLedger::EnergyLedger(std::vector<int> region_bounds)
: region_bounds_(std::move(region_bounds))
{
}

void EnergyLedger::record(const std::string& label, const ParticleContainer& particles) const {
  const int num_regions = std::ssize(region_bounds_) - 1;
  std::vector<double> weights(num_regions, 0.0), energies(num_regions, 0.0);
  for (const Particle& particle : particles) {
    if (not particle.is_alive) continue;
    const auto upper = std::upper_bound(region_bounds_.begin(), region_bounds_.end(), particle.element);
    const int region = std::clamp(static_cast<int>(upper - region_bounds_.begin()) - 1, 0, num_regions - 1);
    weights[region] += particle.weight;
    energies[region] += 0.5 * particle.weight * particle.species.mass * (particle.velocity * particle.velocity);
  }
  std::cout << std::setprecision(12) << "    LEDGER " << label;
  for (int region = 0; region < num_regions; ++region) std::cout << " " << weights[region] << " " << energies[region];
  std::cout << "\n";
}

const EnergyLedger*& activeEnergyLedger() {
  static const EnergyLedger* ledger = nullptr;
  return ledger;
}

void recordEnergyLedger(const std::string& label, const ParticleContainer& particles) {
  if (activeEnergyLedger() != nullptr) activeEnergyLedger()->record(label, particles);
}

} // namespace mfpic
