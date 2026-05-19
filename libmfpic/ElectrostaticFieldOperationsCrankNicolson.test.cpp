#include <libmfpic/Constants.hpp>
#include <libmfpic/DGEulerInitialConditionsFactory.hpp>
#include <libmfpic/DGEulerOperations.hpp>
#include <libmfpic/DGEulerOperationsFactory.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/ElectrostaticFieldOperations.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/LowFidelityState.hpp>
#include <libmfpic/MFEMHelpers.hpp>
#include <libmfpic/Pinning.hpp>
#include <libmfpic/SourcesFactory.hpp>
#include <libmfpic/TestingUtils.hpp>

#include <mfem.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

Species electron_species{.charge = -constants::elementary_charge, .mass = constants::electron_mass};
Species proton_species{.charge = constants::elementary_charge, .mass = constants::proton_mass};

TEST(ElectrostaticFieldOperationsCrankNicolson, test) {
  auto& message = std::cout;
  message << std::setprecision(16);

  constexpr int num_elems = 4;
  constexpr double length = 1.3;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
    num_elems, num_elems, num_elems, mfem::Element::HEXAHEDRON, length, length, length);

  constexpr int dg_order = 0;
  constexpr int num_euler_equations = euler::ConservativeVariables::NUM_VARS;
  Discretization euler_dg_discretization(&mesh, dg_order, FETypes::DG, num_euler_equations);
  constexpr int es_order = 1;
  Discretization es_discretization(&mesh, es_order, FETypes::HGRAD);

  constexpr double number_density_e = 3.e23;
  const mfem::Vector bulk_velocity_e{90.2, 170.1, 315.9};
  constexpr double temperature_e = 101.6;

  constexpr double number_density_p = 8.1e22;
  const mfem::Vector bulk_velocity_p{256.1, 92.5, 165.8};
  constexpr double temperature_p = 83.5;

  const mfem::Vector primitive_e = euler::constructPrimitiveState(number_density_e, bulk_velocity_e, temperature_e);
  const mfem::Vector primitive_p = euler::constructPrimitiveState(number_density_p, bulk_velocity_p, temperature_p);

  const double plasma_frequency_e = sqrt(
    (electron_species.charge * electron_species.charge * number_density_e) / (electron_species.mass * constants::permittivity));

  const double plasma_period_e = 2.0 * M_PI / plasma_frequency_e;

  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters;
  list_of_parameters.push_back(
    std::make_unique<ConstantSourceParameters>(electron_species, number_density_e, temperature_e, bulk_velocity_e));
  list_of_parameters.push_back(
    std::make_unique<ConstantSourceParameters>(proton_species, number_density_p, temperature_p, bulk_velocity_p));

  LowFidelityState low_fidelity_state = buildEulerState(euler_dg_discretization, list_of_parameters);

  const std::vector<Species> species_list = low_fidelity_state.getSpeciesList();
  std::vector<std::vector<std::unique_ptr<DGBC>>> empty_bcs(species_list.size());
  std::vector<std::pair<Species, std::unique_ptr<mfem::VectorCoefficient>>> empty_sources{};
  std::unique_ptr<LowFidelityOperations> dg_euler_operations = buildDGEulerOperations(
    euler_dg_discretization, es_discretization, species_list, empty_bcs, empty_sources);

  ElectrostaticFieldState field_state(es_discretization);

  const double dt = plasma_period_e / 10;
  const auto [updated_low_fidelity_state, updated_field_state] = dg_euler_operations->plasmaOscillate(
    dt, low_fidelity_state, field_state);

  const double total_fluid_energy_initial = dg_euler_operations->computeTotalEnergy(low_fidelity_state);
  const double total_fluid_energy_updated = dg_euler_operations->computeTotalEnergy(updated_low_fidelity_state);
  const double total_kinetic_energy_initial = dg_euler_operations->computeTotalKineticEnergy(low_fidelity_state);
  const double total_kinetic_energy_updated = dg_euler_operations->computeTotalKineticEnergy(updated_low_fidelity_state);
  const double total_internal_energy_initial = total_fluid_energy_initial - total_kinetic_energy_initial;
  const double total_internal_energy_updated = total_fluid_energy_updated - total_kinetic_energy_updated;

  const double relative_tolerance = 1e-14;
  EXPECT_NEAR_RELATIVE(total_internal_energy_initial, total_internal_energy_updated, relative_tolerance);

  auto pinning = std::make_unique<Pinning>();
  ElectrostaticFieldOperations es_field_operations(es_discretization, std::move(pinning));
  const double total_field_energy_initial = es_field_operations.fieldEnergy(field_state);
  const double total_field_energy_updated = es_field_operations.fieldEnergy(updated_field_state);

  const double total_energy_initial = total_fluid_energy_initial + total_field_energy_initial;
  const double total_energy_updated = total_fluid_energy_updated + total_field_energy_updated;

  EXPECT_NEAR_RELATIVE(total_energy_initial, total_energy_updated, relative_tolerance);

  const LowFidelitySpeciesState& updated_species_state_0 = updated_low_fidelity_state.getSpeciesState(0);
  const LowFidelitySpeciesState& updated_species_state_1 = updated_low_fidelity_state.getSpeciesState(1);

  const mfem::GridFunction& updated_species_0 = updated_species_state_0.getGridFunction();
  const mfem::GridFunction& updated_species_1 = updated_species_state_1.getGridFunction();

  const double expected_rho_0 = number_density_e * electron_species.mass;
  const double expected_rho_1 = number_density_p * proton_species.mass;
  // computed in mathematica
  const mfem::Vector expected_p_0({0.00002361634865202284, 0.00003935995377889034, 0.00007301770893578254});
  const mfem::Vector expected_p_1({0.034697314983997654, 0.01253404337795641, 0.022466571774853485});

  const double expected_kinetic_energy_density_0 = euler::kineticEnergyDensity(expected_rho_0, expected_p_0);
  const double expected_kinetic_energy_density_1 = euler::kineticEnergyDensity(expected_rho_1, expected_p_1);
  const double expected_internal_energy_density_0 = euler::getInternalEnergyDensityFromPrimitiveState(primitive_e, electron_species);
  const double expected_internal_energy_density_1 = euler::getInternalEnergyDensityFromPrimitiveState(primitive_p, proton_species);
  const double expected_total_energy_density_0 = expected_internal_energy_density_0 + expected_kinetic_energy_density_0;
  const double expected_total_energy_density_1 = expected_internal_energy_density_1 + expected_kinetic_energy_density_1;

  const mfem::Vector expected_conservative_state_0 = euler::constructConservativeState(
    expected_rho_0, expected_p_0, expected_total_energy_density_0);
  const mfem::Vector expected_conservative_state_1 = euler::constructConservativeState(
    expected_rho_1, expected_p_1, expected_total_energy_density_1);
  mfem::VectorConstantCoefficient expected_coefficient_0(expected_conservative_state_0);
  mfem::VectorConstantCoefficient expected_coefficient_1(expected_conservative_state_1);

  const double error_species_0 = updated_species_0.ComputeL2Error(expected_coefficient_0);
  const double error_species_1 = updated_species_1.ComputeL2Error(expected_coefficient_1);

  constexpr double absolute_tolerance = 1e-16;
  EXPECT_NEAR(0., error_species_0, absolute_tolerance);
  EXPECT_NEAR(0., error_species_1, absolute_tolerance);

  const mfem::GridFunction updated_e_field = updated_field_state.getEFieldGridFunction();
  // computed in mathematica
  const mfem::Vector expected_updated_e_field({2115.150714223121, 14580.40409450579, 27240.2637985167});
  mfem::VectorConstantCoefficient expected_e_field_coefficient(expected_updated_e_field);
  const double error_e_field = updated_e_field.ComputeL2Error(expected_e_field_coefficient);

  constexpr double e_field_tolerance = 1e-10;
  EXPECT_NEAR(0., error_e_field, e_field_tolerance);
}

// CrankNicolson TimeStepper
// plasma oscillation end to end test


}