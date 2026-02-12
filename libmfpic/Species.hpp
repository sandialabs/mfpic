#pragma once

namespace mfpic {

struct Species {
  double charge = 0.0;
  double mass = 1.0;
  double charge_over_mass = charge / mass;
  double specific_heat_ratio = 5. / 3.;

  bool operator==(const Species& other_species) const {
    return
      charge == other_species.charge and
      mass == other_species.mass and
      charge_over_mass == other_species.charge_over_mass and
      specific_heat_ratio == other_species.specific_heat_ratio;
  }
};

} // namespace mfpic
