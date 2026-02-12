#include <libmfpic/DGEulerInitialConditionsFactory.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/LowFidelityState.hpp>
#include <libmfpic/SourcesFactory.hpp>

namespace mfpic {

LowFidelitySpeciesState buildEulerSpeciesState(Discretization& discretization, const SourceParameters& parameters) {
  std::unique_ptr<mfem::VectorCoefficient> initial_conditions_coefficient = parameters.getEulerVectorCoefficient();
  LowFidelitySpeciesState euler_species_state(discretization, parameters.species, *initial_conditions_coefficient);

  return euler_species_state;
}

LowFidelityState buildEulerState(
  Discretization& discretization,
  const std::vector<std::unique_ptr<SourceParameters>>& list_of_parameters)
{
  std::vector<std::pair<Species, std::unique_ptr<mfem::VectorCoefficient>>> species_coefficient_list;
  for (const std::unique_ptr<SourceParameters>& parameters : list_of_parameters) {
    std::unique_ptr<mfem::VectorCoefficient> initial_conditions_coefficient = parameters->getEulerVectorCoefficient();
    species_coefficient_list.emplace_back(
      std::piecewise_construct,
      std::forward_as_tuple(parameters->species),
      std::forward_as_tuple(std::move(initial_conditions_coefficient)));
  }

  LowFidelityState euler_state(discretization, species_coefficient_list);
  return euler_state;
}

}