#pragma once

#include <libmfpic/ParticleContainer.hpp>

#include <string>
#include <vector>

namespace mfpic {

/**
 * @brief Diagnostic that prints, per spatial region of elements, the total macroparticle count weight and kinetic
 *        energy of all live particles, tagged with a label.
 *
 * Output lines have the form "    LEDGER <label> <w_0> <E_0> <w_1> <E_1> ...", one pair per region. Differences between
 * consecutive labels within a timestep attribute energy changes to individual substeps.
 */
class EnergyLedger {
public:
  /// @param[in] region_bounds Element-index boundaries [b_0, b_1, ..., b_R]; region r holds elements in [b_r, b_{r+1}).
  explicit EnergyLedger(std::vector<int> region_bounds);

  void record(const std::string& label, const ParticleContainer& particles) const;

private:
  std::vector<int> region_bounds_;
};

/// Ledger used by the time loop when the deck sets Output: Energy Ledger Regions; null otherwise.
const EnergyLedger*& activeEnergyLedger();

/// Record to the active ledger, if any.
void recordEnergyLedger(const std::string& label, const ParticleContainer& particles);

} // namespace mfpic
