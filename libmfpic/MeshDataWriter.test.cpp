#include <libmfpic/Constants.hpp>
#include <libmfpic/DGEulerInitialConditionsFactory.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/ElectrostaticFieldState.hpp>
#include <libmfpic/LowFidelityState.hpp>
#include <libmfpic/MeshDataWriter.hpp>
#include <libmfpic/MeshFactory.hpp>
#include <libmfpic/SourcesFactory.hpp>
#include <libmfpic/Species.hpp>

#include <gtest/gtest.h>

#include <filesystem>

namespace {

using namespace mfpic;

const Species electron_species{.charge = -constants::elementary_charge, .mass = constants::electron_mass};
const Species proton_species{.charge = constants::elementary_charge, .mass = constants::proton_mass};

TEST(MeshDataWriter, test) {
  const MeshParameters mesh_parameters{
    .mesh_type = "quad",
    .lengths{1., 1.},
    .num_elements{10, 10}
  };
  mfem::Mesh mesh = buildMesh(mesh_parameters);

  const std::string folder_name{"MeshOutput"};
  MeshDataWriter mesh_data_writer(folder_name, mesh);

  constexpr int electrostatic_order = 1;
  Discretization electrostatic_discretization(&mesh, electrostatic_order, FETypes::HGRAD);
  ElectrostaticFieldState electrostatic_field_state(electrostatic_discretization);

  // Nonzero data is being put into ElectrostaticFieldState and LowFidelityState to be manually checked in output but won't
  // be checked specifically in this unit test

  mfem::GridFunction potential(&electrostatic_discretization.getFeSpace());
  mfem::FunctionCoefficient quadratic_function([](const mfem::Vector& x){ return x[0] * x[0] - x[1] * x[1];});
  potential.ProjectCoefficient(quadratic_function);
  electrostatic_field_state.setPotential(potential);

  constexpr int dg_order = 0;
  constexpr int num_equations = 5;
  Discretization dg_discretization(&mesh, dg_order, FETypes::DG, num_equations);

  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters;
  list_of_parameters.push_back(std::make_unique<ConstantSourceParameters>(electron_species, 1e25, 300));
  list_of_parameters.push_back(std::make_unique<ConstantSourceParameters>(proton_species, 2e22, 320));
  std::vector<LowFidelityState> low_fidelity_states = {buildEulerState(dg_discretization, list_of_parameters)};

  mesh_data_writer.output(electrostatic_field_state, low_fidelity_states, 0, 0);
  mesh_data_writer.output(electrostatic_field_state, low_fidelity_states, 1, 1.);

  EXPECT_TRUE(std::filesystem::exists(folder_name));

  const std::string main_file_name = folder_name + "/" + folder_name + ".pvd";
  EXPECT_TRUE(std::filesystem::exists(main_file_name));

  const std::string cycle_0_folder = folder_name + "/Cycle000000";
  const std::string cycle_0_data = cycle_0_folder + "/data.pvtu";
  EXPECT_TRUE(std::filesystem::exists(cycle_0_folder));
  EXPECT_TRUE(std::filesystem::exists(cycle_0_data));

  const std::string cycle_1_folder = folder_name + "/Cycle000001";
  const std::string cycle_1_data = cycle_1_folder + "/data.pvtu";
  EXPECT_TRUE(std::filesystem::exists(cycle_1_folder));
  EXPECT_TRUE(std::filesystem::exists(cycle_1_data));
}

}