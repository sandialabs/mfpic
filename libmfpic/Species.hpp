#pragma once

#include <string>

namespace mfpic {

struct Species {
  double charge = 0.0;
  double mass = 1.0;
  double charge_over_mass = charge / mass;
  double specific_heat_ratio = 5. / 3.;
  std::string name = "";

  auto operator <=>(const Species& other_species) const = default;
};

} // namespace mfpic

namespace std {
  template <typename... Ts>
    static std::size_t hashTuple(const std::tuple<Ts...>& tup) {
        std::size_t seed = 0;

        std::apply([&seed](const auto&... elems) {
            ((seed ^= std::hash<std::decay_t<decltype(elems)>>{}(elems)
                       + 0x9e3779b9 + (seed << 6) + (seed >> 2)), ...);
        }, tup);

        return seed;
  }
  template <>
  struct hash<mfpic::Species> {
    std::size_t operator()(const mfpic::Species& species) const noexcept {
      return hashTuple(std::tie(species.charge, 
                                     species.mass, 
                                     species.charge_over_mass, 
                                     species.specific_heat_ratio, 
                                     species.name));
    }
  };
} // namespace std
