#include <libmfpic/LowFidelityState.hpp>

#include <libmfpic/Discretization.hpp>

namespace mfpic {

LowFidelitySpeciesState::LowFidelitySpeciesState(Discretization& discretization, const Species& species)
  : grid_function_(&discretization.getFeSpace())
  , species_(species)
{
  grid_function_ = 0;
}

LowFidelitySpeciesState::LowFidelitySpeciesState(
  Discretization& discretization,
  const Species& species,
  mfem::VectorCoefficient& coefficient)
  : grid_function_(&discretization.getFeSpace())
  , species_(species)
{
  grid_function_.ProjectCoefficient(coefficient);
}

LowFidelityState::LowFidelityState(Discretization& discretization, const std::vector<Species>& species_list) {
  for (const Species& species : species_list) {
    species_states_.emplace_back(discretization, species);
  }
}

LowFidelityState::LowFidelityState(
  Discretization& discretization,
  const std::vector<std::pair<Species, std::unique_ptr<mfem::VectorCoefficient>>>& species_coefficient_list)
{
  for (const auto& [species, coefficient] : species_coefficient_list) {
    species_states_.emplace_back(discretization, species, *coefficient);
  }
}

std::vector<Species> LowFidelityState::getSpeciesList() const {
  std::vector<Species> species_list;
  for (const LowFidelitySpeciesState& species_state : species_states_) {
    species_list.push_back(species_state.getSpecies());
  }
  return species_list;
}

void LowFidelityState::addScaledState(const double scale, const LowFidelityState& state_to_sum) {
  assert(numSpecies() == state.numSpecies());

  for (int i_species = 0; i_species < numSpecies(); ++i_species) {
    LowFidelitySpeciesState& species_state = getSpeciesState(i_species);
    const LowFidelitySpeciesState& species_state_to_sum = state_to_sum.getSpeciesState(i_species);

    mfem::GridFunction& species_grid_function = species_state.getGridFunction();
    const mfem::GridFunction& species_grid_function_to_sum = species_state_to_sum.getGridFunction();

    assert(species_grid_function.Size() == species_grid_function_to_sum.Size());

    species_grid_function.Add(scale, species_grid_function_to_sum);
  }
}

}