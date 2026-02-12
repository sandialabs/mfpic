#pragma once

#include <libmfpic/Species.hpp>
#include <mfem.hpp>

namespace mfpic {

class Discretization;
/**
 * @brief class to store all information needed to define the low fidelity state for a single species
 */
class LowFidelitySpeciesState {
public:
  /**
   * @brief Construct a new LowFidelitySpeciesState object initialized to zero
   * 
   * @param discretization - the discretization of the species, this must be non const because MFEM grid functions requires a
   *  non const finite element space
   * @param species - the data defining the species
   */
  LowFidelitySpeciesState(Discretization& discretization, const Species& species);

  /**
   * @brief Construct a LowFidelityState initialized to a vector coefficient
   *
   * @param discretization - the discretization of the species, this must be non const because MFEM grid functions requires a
   *  non const finite element space
   * @param species - the data defining the species
   * @param coefficient - coefficient to project onto grid function, this must be non const because MFEM grid functions can only
   *  project non const coefficients.
   */
  LowFidelitySpeciesState(
    Discretization& discretization,
    const Species& species,
    mfem::VectorCoefficient& coefficient);

  Species getSpecies() const { return species_; }

  mfem::GridFunction& getGridFunction() { return grid_function_; }
  const mfem::GridFunction& getGridFunction() const { return grid_function_; }

private:
  mfem::GridFunction grid_function_;

  const Species species_;
};

class LowFidelityState {

public:
  /**
   * @brief Construct a LowFidelityState initialized to zero
   *
   * @param discretization - discretization of each species
   * @param list_of_species - List of \ref Species
   */
  LowFidelityState(Discretization& discretization, const std::vector<Species>& list_of_species);

  /**
   * @brief Construct a LowFidelityState initialized to a list of coefficients
   *
   * @param discretization - discretization of each species
   * @param species_coefficient_list - list of species and the coefficient defining the initial condition for that species
   */
  LowFidelityState(
    Discretization& discretization,
    const std::vector<std::pair<Species, std::unique_ptr<mfem::VectorCoefficient>>>& species_coefficient_list);

  /// Number of species in the state
  int numSpecies() const { return std::ssize(species_states_); };
  std::vector<Species> getSpeciesList();

  LowFidelitySpeciesState& getSpeciesState(const int i_species) { return species_states_[i_species]; }
  const LowFidelitySpeciesState& getSpeciesState(const int i_species) const { return species_states_[i_species]; }

private:
  std::vector<LowFidelitySpeciesState> species_states_;
};

} // namespace mfpic
