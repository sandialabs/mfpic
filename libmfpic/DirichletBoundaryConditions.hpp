#pragma once

#include <mfem/mfem.hpp>

namespace mfpic {

class DirichletBoundaryConditions {
public:
  virtual ~DirichletBoundaryConditions() = default;

  virtual mfem::Array<int> getDirichletBoundaryDofIndices() const = 0;

  virtual void applyBoundaryConditions(mfem::GridFunction& solution) const = 0;
};

}