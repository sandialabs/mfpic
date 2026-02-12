#include <libmfpic/Constants.hpp>
#include <libmfpic/DGEulerInitialConditionsFactory.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/LowFidelityState.hpp>
#include <libmfpic/SourcesFactory.hpp>
#include <libmfpic/Species.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

Species electron_species{.charge = -constants::elementary_charge, .mass = constants::electron_mass};
Species proton_species{.charge = constants::elementary_charge, .mass = constants::proton_mass};

TEST(DGEulerInitialConditionsFactory, buildEulerSpeciesState) {
  constexpr int num_elems = 10;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);

  constexpr int dg_order = 0;
  constexpr int num_equations = euler::ConservativeVariables::NUM_VARS;
  Discretization dg_discretization(&mesh, dg_order, FETypes::DG, num_equations);

  constexpr double number_density = 1.e16;
  const mfem::Vector bulk_velocity{12.3, 45.6, 78.9};
  constexpr double temperature = 300;

  ConstantSourceParameters parameters(electron_species, number_density, temperature, bulk_velocity);

  LowFidelitySpeciesState euler_species_state = buildEulerSpeciesState(dg_discretization, parameters);

  mfem::GridFunction species_grid_function = euler_species_state.getGridFunction();
  EXPECT_EQ(num_elems * num_equations, species_grid_function.Size());

  mfem::Vector primitive_state = euler::constructPrimitiveState(number_density, bulk_velocity, temperature);
  mfem::Vector conservative_state = euler::convertFromPrimitiveToConservative(primitive_state, electron_species);

  mfem::VectorConstantCoefficient exact_coefficient(conservative_state);
  const double l2_error = species_grid_function.ComputeL2Error(exact_coefficient);

  EXPECT_DOUBLE_EQ(0., l2_error);
}

TEST(DGEulerInitialConditionsFactory, buildEulerState) {
  constexpr int num_elems = 15;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);

  constexpr int dg_order = 0;
  constexpr int num_equations = euler::ConservativeVariables::NUM_VARS;
  Discretization dg_discretization(&mesh, dg_order, FETypes::DG, num_equations);

  constexpr double number_density_0 = 2.e17;
  const mfem::Vector bulk_velocity_0{21.9, 87.6, 54.3};
  constexpr double temperature_0 = 287;

  constexpr double number_density_1 = 9.e16;
  const mfem::Vector bulk_velocity_1{35.7, 46.8, 57.9};
  constexpr double temperature_1 = 291;

  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters;
  list_of_parameters.push_back(
    std::make_unique<ConstantSourceParameters>(electron_species, number_density_0, temperature_0, bulk_velocity_0));
  list_of_parameters.push_back(
    std::make_unique<ConstantSourceParameters>(proton_species, number_density_1, temperature_1, bulk_velocity_1));

  LowFidelityState euler_state = buildEulerState(dg_discretization, list_of_parameters);

  EXPECT_EQ(2, euler_state.numSpecies());

  const LowFidelitySpeciesState& euler_species_state_0 = euler_state.getSpeciesState(0);
  EXPECT_EQ(electron_species, euler_species_state_0.getSpecies());

  const mfem::GridFunction& species_grid_function_0 = euler_species_state_0.getGridFunction();

  mfem::Vector primitive_state_0 = euler::constructPrimitiveState(number_density_0, bulk_velocity_0, temperature_0);
  mfem::Vector conservative_state_0 = euler::convertFromPrimitiveToConservative(primitive_state_0, electron_species);
  mfem::VectorConstantCoefficient exact_coefficient_0(conservative_state_0);

  const double l2_error_0 = species_grid_function_0.ComputeL2Error(exact_coefficient_0);
  EXPECT_DOUBLE_EQ(0., l2_error_0);

  const LowFidelitySpeciesState& euler_species_state_1 = euler_state.getSpeciesState(1);
  EXPECT_EQ(proton_species, euler_species_state_1.getSpecies());

  const mfem::GridFunction& species_grid_function_1 = euler_species_state_1.getGridFunction();

  mfem::Vector primitive_state_1 = euler::constructPrimitiveState(number_density_1, bulk_velocity_1, temperature_1);
  mfem::Vector conservative_state_1 = euler::convertFromPrimitiveToConservative(primitive_state_1, proton_species);
  mfem::VectorConstantCoefficient exact_coefficient_1(conservative_state_1);

  const double l2_error_1 = species_grid_function_1.ComputeL2Error(exact_coefficient_1);
  EXPECT_DOUBLE_EQ(0., l2_error_1);
}

}