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

std::vector<Species> LowFidelityState::getSpeciesList() {
  std::vector<Species> species_list;
  for (const LowFidelitySpeciesState& species_state : species_states_) {
    species_list.push_back(species_state.getSpecies());
  }
  return species_list;
}
}