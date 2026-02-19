#include <libmfpic/DGEulerBoundaryConditions.hpp>
#include <libmfpic/DGEulerInitialConditionsFactory.hpp>
#include <libmfpic/DGGhostBoundaryIntegrator.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/DGEulerOperationsFactory.hpp>
#include <libmfpic/LowFidelityOperations.hpp>
#include <libmfpic/LowFidelityState.hpp>
#include <libmfpic/SourcesFactory.hpp>
#include <libmfpic/Species.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

constexpr Species species_0{.charge =  1.0, .mass = 2.0, .specific_heat_ratio = 1.4};
constexpr Species species_1{.charge = -1.0, .mass = 3.0, .specific_heat_ratio = 2.8};
constexpr Species species_2{.charge =  4.0, .mass = 6.0, .specific_heat_ratio = 4.2};
constexpr Species species_3{.charge = -2.0, .mass = 1.5, .specific_heat_ratio = 5.6};

const std::vector<Species> species_list{species_0,species_1,species_2,species_3};

const std::vector<std::unique_ptr<DGGhostBC>> empty_bcs{};

TEST(DGEulerOperationsFactory, ErrorsOutIfDiscretizationIsWrong) {
  constexpr int nx = 5;
  constexpr mfem::real_t length = 1.0;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(nx, nx, nx, mfem::Element::HEXAHEDRON, length, length, length);

  constexpr int basis_order = 1;
  int num_equations = 5;
  Discretization charge_discretization(&mesh, basis_order, FETypes::HGRAD);
  Discretization discretization(&mesh, basis_order, FETypes::HGRAD, num_equations);

  EXPECT_DEATH(buildDGEulerOperations(discretization, charge_discretization, species_list, empty_bcs),
               "Fluid discretization must be DG.");
}

TEST(DGEulerOperationsFactory, BasicChecksOnOperationsAndState) {
  constexpr int nx = 5;
  constexpr mfem::real_t length = 1.0;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(nx, nx, nx, mfem::Element::HEXAHEDRON, length, length, length);

  constexpr int basis_order = 1;
  int num_equations = 5;
  Discretization charge_discretization(&mesh, basis_order, FETypes::HGRAD);
  Discretization discretization(&mesh, basis_order, FETypes::DG, num_equations);

  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters;
  for (const Species& species : species_list) {
    const double number_density = 1000;
    const double temperature = 300;

    list_of_parameters.push_back(std::make_unique<ConstantSourceParameters>(species, number_density, temperature));
  }
  LowFidelityState dg_euler_state = buildEulerState(discretization, list_of_parameters);

  constexpr int boundary_attribute = 3;
  auto bc = std::make_unique<DGEulerReflectingBC>(boundary_attribute, mesh);
  std::vector<std::unique_ptr<DGGhostBC>> bcs;
  bcs.push_back(std::move(bc));

  std::unique_ptr<LowFidelityOperations> dg_euler_operations = buildDGEulerOperations(
    discretization, charge_discretization, species_list, bcs);

  EXPECT_NO_FATAL_FAILURE(dg_euler_operations->move(1.0, dg_euler_state));
}

} // namespace
