#pragma once

#include <libmfpic/Species.hpp>

#include <mfem/mfem.hpp>

namespace mfpic {

struct Particle {
  mfem::Vector position = mfem::Vector({0.0, 0.0, 0.0});
  mfem::Vector velocity = mfem::Vector({0.0, 0.0, 0.0});
  int element = 0;
  Species species = Species{};
  double weight = 0.0;
  bool is_alive = true;
};

} // namespace mfpic
