#include <libmfpic/Constants.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/LowFidelityState.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

const Species electron_species{.charge = -constants::elementary_charge, .mass = constants::electron_mass};
const Species proton_species{.charge = constants::elementary_charge, .mass = constants::proton_mass};

TEST(LowFidelityState, constructLowFidelitySpeciesState) {
  constexpr int num_elems = 10;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);

  constexpr int dg_order = 0;
  constexpr int num_equations = euler::ConservativeVariables::NUM_VARS;
  Discretization dg_discretization(&mesh, dg_order, FETypes::DG, num_equations);

  const LowFidelitySpeciesState species_state(dg_discretization, electron_species);

  EXPECT_EQ(electron_species, species_state.getSpecies());

  const mfem::GridFunction& species_grid_function = species_state.getGridFunction();
  EXPECT_EQ(&dg_discretization.getFeSpace(), species_grid_function.FESpace());

  const double norm = species_grid_function.Norml2();
  EXPECT_EQ(0., norm);
}

TEST(LowFidelityState, constructLowFidelitySpeciesStateWithVectorCoefficient) {
  constexpr int num_elems = 10;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);

  constexpr int dg_order = 0;
  constexpr int num_equations = 1;
  Discretization dg_discretization(&mesh, dg_order, FETypes::DG, num_equations);

  auto function = [&](const mfem::Vector&, mfem::Vector& y) {
    y[0] = 2.3;
  };
  mfem::VectorFunctionCoefficient vector_coefficient(num_equations, function);

  const LowFidelitySpeciesState species_state(dg_discretization, electron_species, vector_coefficient);

  EXPECT_EQ(electron_species, species_state.getSpecies());

  const mfem::GridFunction& species_grid_function = species_state.getGridFunction();
  EXPECT_EQ(&dg_discretization.getFeSpace(), species_grid_function.FESpace());

  const double error = species_grid_function.ComputeL2Error(vector_coefficient);
  EXPECT_DOUBLE_EQ(0., error);
}

TEST(LowFidelityState, copyConstructorAllocatesNewMemoryLowFidelitySpeciesState) {
  constexpr int num_elems = 10;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);

  constexpr int dg_order = 0;
  constexpr int num_equations = euler::ConservativeVariables::NUM_VARS;
  Discretization dg_discretization(&mesh, dg_order, FETypes::DG, num_equations);

  LowFidelitySpeciesState species_state(dg_discretization, electron_species);
  mfem::GridFunction& grid_function_original = species_state.getGridFunction();
  grid_function_original.Randomize();

  LowFidelitySpeciesState species_state_copy(species_state);
  const mfem::GridFunction& grid_function_copy = species_state_copy.getGridFunction();

  EXPECT_NE(grid_function_original.GetData(), grid_function_copy.GetData());
  for (int i = 0; i < grid_function_original.Size(); ++i) {
    EXPECT_EQ(grid_function_original[i], grid_function_copy[i]);
  }
}

TEST(LowFidelityState, constructLowFidelityState) {
  constexpr int num_elems = 10;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);

  constexpr int dg_order = 0;
  constexpr int num_equations = euler::ConservativeVariables::NUM_VARS;
  Discretization dg_discretization(&mesh, dg_order, FETypes::DG, num_equations);

  const std::vector<Species> species_list{electron_species, proton_species};

  const LowFidelityState state(dg_discretization, species_list);

  EXPECT_EQ(2, state.numSpecies());
  for (int i_species = 0; i_species < 2; ++i_species) {
    const LowFidelitySpeciesState& species_state = state.getSpeciesState(i_species);
    EXPECT_EQ(species_list[i_species], species_state.getSpecies());

    const mfem::GridFunction& grid_function = species_state.getGridFunction();
    const double norm = grid_function.Norml2();
    EXPECT_EQ(0., norm);
  }
}

TEST(LowFidelityState, constructLowFidelityStateWithVectorCoefficients) {
  constexpr int num_elems = 10;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);

  constexpr int dg_order = 0;
  constexpr int num_equations = 1;
  Discretization dg_discretization(&mesh, dg_order, FETypes::DG, num_equations);

  auto function_0 = [&](const mfem::Vector&, mfem::Vector& y) { y[0] = 2.3; };
  auto function_1 = [&](const mfem::Vector&, mfem::Vector& y) { y[0] = 3.4; };

  auto vector_coefficient_0 = std::make_unique<mfem::VectorFunctionCoefficient>(num_equations, function_0);
  auto vector_coefficient_1 = std::make_unique<mfem::VectorFunctionCoefficient>(num_equations, function_1);

  std::vector<std::pair<Species, std::unique_ptr<mfem::VectorCoefficient>>> species_coefficient_list;
  species_coefficient_list.push_back(std::make_pair(electron_species, std::move(vector_coefficient_0)));
  species_coefficient_list.push_back(std::make_pair(proton_species, std::move(vector_coefficient_1)));

  const LowFidelityState state(dg_discretization, species_coefficient_list);

  EXPECT_EQ(2, state.numSpecies());
  for (int i_species = 0; i_species < 2; ++i_species) {
    const LowFidelitySpeciesState& species_state = state.getSpeciesState(i_species);

    const auto& [species, vector_coefficient] = species_coefficient_list[i_species];

    EXPECT_EQ(species, species_state.getSpecies());

    const mfem::GridFunction& grid_function = species_state.getGridFunction();
    const double error = grid_function.ComputeL2Error(*vector_coefficient);
    EXPECT_EQ(0., error);
  }
}

TEST(LowFidelityState, copyConstructorAllocatesNewMemoryLowFidelityState) {
  constexpr int num_elems = 10;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);

  constexpr int dg_order = 0;
  constexpr int num_equations = euler::ConservativeVariables::NUM_VARS;
  Discretization dg_discretization(&mesh, dg_order, FETypes::DG, num_equations);

  const Species species{.charge = -constants::elementary_charge, .mass = constants::electron_mass};
  const std::vector<Species> species_list{species};

  LowFidelityState state_orginal(dg_discretization, species_list);
  LowFidelitySpeciesState& species_state_original = state_orginal.getSpeciesState(0);
  mfem::GridFunction& species_grid_function_original = species_state_original.getGridFunction();
  species_grid_function_original.Randomize();

  LowFidelityState state_copy(state_orginal);
  LowFidelitySpeciesState& species_state_copy = state_copy.getSpeciesState(0);
  mfem::GridFunction& species_grid_function_copy = species_state_copy.getGridFunction();

  EXPECT_NE(&species_state_original, &species_state_copy);
  EXPECT_NE(&species_grid_function_original, &species_grid_function_copy);
  EXPECT_NE(species_grid_function_original.GetData(), species_grid_function_copy.GetData());

  for (int i = 0; i < species_grid_function_original.Size(); ++i) {
    EXPECT_EQ(species_grid_function_original[i], species_grid_function_copy[i]);
  }
}

}