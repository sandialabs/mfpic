#include <libmfpic/TextDataWriter.hpp>

#include <libmfpic/DGBC.hpp>
#include <libmfpic/DGEulerInitialConditionsFactory.hpp>
#include <libmfpic/DGEulerOperationsFactory.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/ElectrostaticFieldOperations.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/LowFidelityOperations.hpp>
#include <libmfpic/LowFidelityState.hpp>
#include <libmfpic/MeshFactory.hpp>
#include <libmfpic/MeshUtilities.hpp>
#include <libmfpic/Pinning.hpp>
#include <libmfpic/SourcesFactory.hpp>

#include <gtest/gtest.h>

#include <filesystem>
#include <ranges>

namespace {

using namespace mfpic;

std::vector<std::pair<Species, std::unique_ptr<mfem::VectorCoefficient>>> empty_sources{};

TEST(TextDataWriter, constructorAndOutput) {
  constexpr int num_low_fidelity_models = 2;

  const MeshParameters mesh_parameters{
    .mesh_type = "quad",
    .lengths{1., 1.},
    .num_elements{10, 10}
  };
  mfem::Mesh mesh = buildMesh(mesh_parameters);

  constexpr int electrostatic_order = 1;
  Discretization electrostatic_discretization(&mesh, electrostatic_order, FETypes::HGRAD);
  constexpr int dg_order = 0;
  Discretization dg_discretization(&mesh, dg_order, FETypes::DG, euler::ConservativeVariables::NUM_VARS);

  auto pinning = std::make_unique<Pinning>();
  ElectrostaticFieldOperations electrostatic_field_operations(electrostatic_discretization, std::move(pinning));
  ElectrostaticFieldState particle_electrostatic_field_state(electrostatic_discretization);

  constexpr int num_species = 1;
  const std::vector<Species> species_list(num_species);
  std::vector<std::vector<std::unique_ptr<DGBC>>> empty_bcs(species_list.size());

  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters;
  constexpr double number_density = 1.;
  constexpr double temperature = 1.;
  list_of_parameters.push_back(std::make_unique<ConstantSourceParameters>(species_list[0], number_density, temperature));

  std::vector<ElectrostaticFieldState> low_fidelity_field_states;
  std::vector<LowFidelityState> low_fidelity_states;
  std::vector<std::unique_ptr<LowFidelityOperations>> low_fidelity_operations;

  for (int i_low_fidelity_model = 0; i_low_fidelity_model < num_low_fidelity_models; ++i_low_fidelity_model) {
    low_fidelity_field_states.emplace_back(electrostatic_discretization);

    LowFidelityState dg_euler_state = buildEulerState(dg_discretization, list_of_parameters);
    low_fidelity_states.push_back(dg_euler_state);

    std::unique_ptr<LowFidelityOperations> dg_euler_operations = buildDGEulerOperations(
      dg_discretization,
      electrostatic_discretization,
      species_list,
      empty_bcs,
      empty_sources);
    low_fidelity_operations.push_back(std::move(dg_euler_operations));
  }

  TextDataWriter text_data_writer(num_low_fidelity_models);

  constexpr int i_timestep = 2;
  constexpr double time = 1.3;
  constexpr double timestep_size = 0.65;
  const double smallest_cell_lengthscale = getSmallestCellLengthscale(mesh);
  text_data_writer.output(
    particle_electrostatic_field_state,
    low_fidelity_field_states,
    electrostatic_field_operations,
    low_fidelity_states,
    low_fidelity_operations,
    timestep_size,
    smallest_cell_lengthscale,
    i_timestep,
    time);

  std::vector<std::string> filenames{"output.csv"};
  for (int i_low_fidelity_model = 0; i_low_fidelity_model < num_low_fidelity_models; ++i_low_fidelity_model) {
    const std::string filename = "output_lf_" + std::to_string(i_low_fidelity_model) + ".csv";
    filenames.push_back(filename);
  }

  for (const std::string& filename : filenames) {
    EXPECT_TRUE(std::filesystem::exists(filename));
    std::ifstream csv_file(filename);

    std::string line_0;
    EXPECT_TRUE(std::getline(csv_file, line_0));
    EXPECT_EQ('#', line_0[0]);
    // remove comment and first space
    line_0.erase(0, 2);

    std::string line_1;
    EXPECT_TRUE(std::getline(csv_file, line_1));

    std::unordered_map<std::string, double> header_to_value;
    auto line_0_split = line_0 | std::views::split(' ');
    auto line_1_split = line_1 | std::views::split(' ');
    for (auto [header, value] : std::views::zip(line_0_split, line_1_split)) {
      std::string header_as_string{std::string_view(header)};
      std::string value_as_string{std::string_view(value)};
      header_to_value[header_as_string] = std::stod(value_as_string);
    }

    EXPECT_EQ(i_timestep, header_to_value.at("Time_Step"));
    EXPECT_EQ(time, header_to_value.at("Time"));
    constexpr double expected_field_energy = 0.;
    EXPECT_EQ(expected_field_energy, header_to_value.at("Field_Energy"));
    // It isn't possible to set up a fluid with zero energy, so instead create fluid with very small energy and check that output
    // value is small and nonnegative.
    constexpr double max_fluid_energy = 1e-20;
    constexpr double min_fluid_energy = 0.;
    EXPECT_TRUE(header_to_value.at("Total_Fluid_Energy") < max_fluid_energy);
    EXPECT_TRUE(header_to_value.at("Total_Fluid_Energy") >= min_fluid_energy);
    constexpr double expected_total_fluid_kinetic_energy = 0.;
    EXPECT_EQ(expected_total_fluid_kinetic_energy, header_to_value.at("Total_Fluid_Kinetic_Energy"));

    std::string line_2;
    EXPECT_FALSE(std::getline(csv_file, line_2));
  }
}

}
