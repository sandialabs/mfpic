#pragma once

#include <mfem/mfem.hpp>

namespace mfpic {

class ElectromagneticFieldsEvaluator {
public:
  virtual mfem::Vector getEFieldAt(const mfem::Vector& position, const int element) const = 0;

  virtual mfem::Vector getBFieldAt(const mfem::Vector& position, const int element) const = 0;
};

} // namespace mfpic
