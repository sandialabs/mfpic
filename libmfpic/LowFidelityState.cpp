#include <libmfpic/LowFidelityState.hpp>

#include <libmfpic/Discretization.hpp>

#include <ranges>

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

LowFidelityState::iterator LowFidelityState::begin() { return species_states_.begin(); }
LowFidelityState::const_iterator LowFidelityState::begin() const { return species_states_.begin(); }
LowFidelityState::iterator LowFidelityState::end() { return species_states_.end(); }
LowFidelityState::const_iterator LowFidelityState::end() const { return species_states_.end(); }

int LowFidelityState::numSpecies() const { return std::ssize(species_states_); }

std::vector<Species> LowFidelityState::getSpeciesList() const {
  std::vector<Species> species_list;
  for (const LowFidelitySpeciesState& species_state : species_states_) {
    species_list.push_back(species_state.getSpecies());
  }
  return species_list;
}

LowFidelitySpeciesState& LowFidelityState::getSpeciesState(const int i_species) { return species_states_[i_species]; }
const LowFidelitySpeciesState& LowFidelityState::getSpeciesState(const int i_species) const { return species_states_[i_species]; }

void LowFidelityState::addScaledState(const double scale, const LowFidelityState& state_to_sum) {
  assert(numSpecies() == state.numSpecies());

  for (auto&& [species_state, species_state_to_sum] : std::views::zip(*this, state_to_sum)) {
    mfem::GridFunction& species_grid_function = species_state.getGridFunction();
    const mfem::GridFunction& species_grid_function_to_sum = species_state_to_sum.getGridFunction();

    assert(species_grid_function.Size() == species_grid_function_to_sum.Size());

    species_grid_function.Add(scale, species_grid_function_to_sum);
  }
}

}