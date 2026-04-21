#pragma once

#include <memory>
#include <mfem.hpp>

namespace mfpic {

struct DGBC;
class Discretization;
class LowFidelityOperations;
class LowFidelityState;
struct Species;

/**
  * @brief Build DGEulerOperations
  *
  * @param dg_discretization \ref Discretization object for the fluids.
  * @param charge_discretization \ref Discretization object for the integrated charge
  * @param list_of_species Species to be evolved with this model
  * @param bcs Boundary conditions
  *
  * @note We currently require that each species use the same \p dg_discretization and
  * the same boundary conditions
  */

std::unique_ptr<LowFidelityOperations>
buildDGEulerOperations(
  Discretization & dg_discretization,
  Discretization & charge_discretization,
  const std::vector<Species>& species_list,
  std::vector<std::vector<std::unique_ptr<DGBC>>> & bcs,
  std::vector<std::pair<Species, std::unique_ptr<mfem::VectorCoefficient>>> & sources)
  ;

} // namespace
