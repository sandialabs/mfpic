#pragma once

#include <memory>
#include <mfem.hpp>

namespace mfpic {

class DGGhostBC;
class Discretization;
class LowFidelityOperations;
class LowFidelityState;
class Species;

/**
  * @brief Build DGEulerOperations
  *
  * @param dg_discretization \ref Discretization object for the fluids.
  * @param charge_discretization \ref Discretization object for the integrated charge
  * @param list_of_species Species to be evolved with this model
  * @param ghost_bcs Boundary conditions applied with the ghost state approach
  *
  * @note We currently require that each species use the same \p dg_discretization and
  * the same boundary conditions
  */

std::unique_ptr<LowFidelityOperations>
buildDGEulerOperations(
  Discretization & dg_discretization,
  Discretization & charge_discretization,
  const std::vector<Species>& species_list,
  const std::vector<std::unique_ptr<DGGhostBC>> & ghost_bcs);

} // namespace
