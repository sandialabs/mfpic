#pragma once

#include <mfem/mfem.hpp>

#include <random>

namespace mfpic {

template <std::uniform_random_bit_generator Generator>
mfem::Vector generateMaxwellianVelocity(
  mfem::Vector bulk_velocity,
  double temperature,
  Generator& generator
);

} // namespace mfpic

#include <libmfpic/GenerateMaxwellianVelocity.tpp>
