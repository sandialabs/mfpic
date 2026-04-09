#include <libmfpic/DGBC.hpp>
#include <libmfpic/DGEulerAssembly.hpp>
#include <libmfpic/DGEulerOperationsFactory.hpp>
#include <libmfpic/DGEulerOperations.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/Errors.hpp>
#include <libmfpic/LowFidelityOperations.hpp>
#include <libmfpic/LowFidelityState.hpp>
#include <memory>
#include <mfem/fem/hyperbolic.hpp>
#include <ranges>

namespace mfpic {

std::unique_ptr<LowFidelityOperations> buildDGEulerOperations(
  Discretization & dg_discretization,
  Discretization & charge_discretization,
  const std::vector<Species>& species_list,
  std::vector<std::vector<std::unique_ptr<DGBC>>> & bcs)
{
  if (dg_discretization.getElementType() != FETypes::DG) {
    std::string error_message = "Fluid discretization must be DG.";
    errorWithDeveloperMessage(error_message);
  }

  std::vector<std::shared_ptr<DGEulerAssembly>> dg_assemblers;
  for (const auto & [i_species, species] : std::views::enumerate(species_list)) {
    dg_assemblers.push_back(std::make_shared<DGEulerAssembly>(dg_discretization.getFeSpace(), species));
    auto & species_bcs = bcs[i_species];
    for (size_t ibc = 0; ibc < species_bcs.size(); ++ibc) {
      dg_assemblers.back()->addBoundaryCondition(std::move(species_bcs[ibc]));
    }
  }

  auto dg_operations = std::make_unique<DGEulerOperations>(charge_discretization, dg_assemblers);
  return dg_operations;
}


} // namespace
