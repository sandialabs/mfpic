#pragma once

#include <memory>
#include <vector>

namespace mfpic {

class Discretization;
struct SourceParameters;
class LowFidelitySpeciesState;
class LowFidelityState;

LowFidelitySpeciesState buildEulerSpeciesState(Discretization& discretization, const SourceParameters& parameters);

/**
 * @brief build an Euler state for one or more species
 * 
 * @param discretization - the discretization used for this Euler species
 * @param list_of_parameters - parameters defining the Euler fluid for each species
 * @return LowFidelityState - the euler state that was built
 */
LowFidelityState buildEulerState(
  Discretization& discretization,
  const std::vector<std::unique_ptr<SourceParameters>>& list_of_parameters);

}