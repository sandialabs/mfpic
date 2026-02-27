#include "libmfpic/IntegratedCharge.hpp"
#include <libmfpic/Constants.hpp>
#include <libmfpic/DirichletBoundaryConditionsConstant.hpp>
#include <libmfpic/ElectrostaticFieldOperationsWithBoltzmannElectrons.hpp>
#include <libmfpic/ElectrostaticFieldState.hpp>
#include <libmfpic/MeshUtilities.hpp>
#include <libmfpic/Pinning.hpp>

#include <gtest/gtest.h>
#include <mfem/mfem.hpp>

using namespace mfpic;

TEST(ElectrostaticFieldOperationsWithBoltzmannElectrons, ZeroChargeAndZeroReferenceDensityWithNaturalBoundariesGivesZeroPotential) {
  constexpr int num_elems = 10;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);
  constexpr int hgrad_order = 1;
  Discretization es_discretization(&mesh, hgrad_order);
  constexpr double reference_number_density = 0.0;
  constexpr double temperature = 1.0;
  mfpic::ElectrostaticFieldOperationsWithBoltzmannElectrons es_field_operations(
    es_discretization,
    std::make_unique<Pinning>(),
    reference_number_density,
    temperature
  );

  ElectrostaticFieldState es_field_state(es_discretization);
  IntegratedCharge integrated_charge(es_discretization);
  integrated_charge.setIntegratedChargeValue(0.);
  es_field_operations.fieldSolve(es_field_state, integrated_charge);

  mfem::GridFunction potential = es_field_state.getPotential();
  mfem::ConstantCoefficient exact_potential(0.);
  const double l2_error = potential.ComputeL2Error(exact_potential);
  constexpr double expected_error = 0.;
  EXPECT_DOUBLE_EQ(expected_error, l2_error);
}

TEST(ElectrostaticFieldOperationsWithBoltzmannElectrons, BoltzmannElectronsSolverGivesSameAnswerAsStandardSolverWithZeroReferenceDensity) {
  constexpr int num_elems = 10;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);
  constexpr int hgrad_order = 1;
  Discretization es_discretization(&mesh, hgrad_order);
  constexpr double reference_number_density = 0.0;
  constexpr double temperature = 100.0;

  IntegratedCharge integrated_charge(es_discretization);
  constexpr double integrated_charge_value = constants::permittivity;
  integrated_charge.setIntegratedChargeValue(integrated_charge_value);

  mfpic::ElectrostaticFieldOperations operations_without_boltzmann_electrons(
    es_discretization,
    std::make_unique<Pinning>()
  );
  ElectrostaticFieldState state_without_boltzmann_electrons(es_discretization);
  operations_without_boltzmann_electrons.fieldSolve(state_without_boltzmann_electrons, integrated_charge);

  mfpic::ElectrostaticFieldOperationsWithBoltzmannElectrons operations_with_boltzmann_electrons(
    es_discretization,
    std::make_unique<Pinning>(),
    reference_number_density,
    temperature
  );
  ElectrostaticFieldState state_with_boltzmann_electrons(es_discretization);
  operations_with_boltzmann_electrons.fieldSolve(state_with_boltzmann_electrons, integrated_charge);

  const mfem::GridFunction potential_without_boltzmann_electrons = state_without_boltzmann_electrons.getPotential();
  const mfem::GridFunction potential_with_boltzmann_electrons = state_with_boltzmann_electrons.getPotential();
  EXPECT_DOUBLE_EQ(potential_without_boltzmann_electrons.DistanceTo(potential_with_boltzmann_electrons), 0.0);
}

// Mathematica
//
// solution :=
//  NDSolve[{-u''[x] + Exp[u[x]] == 0, u[0] == 0, u[1] == 0}, u[x], x]
// nodes := Subdivide[10]
// (u[x] /. solution)[[1]] /. x -> nodes
//
// {0., -0.0414356, -0.0732684, -0.0957998, -0.109238, -0.113704,
// -0.109238, -0.0957998, -0.0732684, -0.0414356, 3.27626*10^-9}
TEST(ElectrostaticFieldOperationsWithBoltzmannElectrons, ZeroChargeWithDirichletBoundaries) {
  constexpr int num_elems = 10;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);
  constexpr int hgrad_order = 1;
  Discretization es_discretization(&mesh, hgrad_order);
  std::unordered_map<std::string, int> side_name_to_boundary_attribute = getSideNameToBoundaryAttributeForInlineMeshes(
    mesh.Dimension());
  const int left_boundary_attribute = side_name_to_boundary_attribute["left"];
  const int right_boundary_attribute = side_name_to_boundary_attribute["right"];
  std::unordered_map<int, double> boundary_attribute_to_dirichlet_value{
    {left_boundary_attribute, 0.},
    {right_boundary_attribute, 0.}};
  auto dirichlet_bcs = std::make_unique<DirichletBoundaryConditionsConstant>(
    boundary_attribute_to_dirichlet_value, es_discretization);
  constexpr double reference_number_density = constants::permittivity / constants::elementary_charge;
  constexpr double temperature = constants::elementary_charge / constants::boltzmann_constant;
  ElectrostaticFieldOperationsWithBoltzmannElectrons es_field_operations(
    es_discretization,
    std::move(dirichlet_bcs),
    reference_number_density,
    temperature
  );
  ElectrostaticFieldState es_field_state(es_discretization);

  IntegratedCharge integrated_charge(es_discretization);
  integrated_charge.setIntegratedChargeValue(0.);

  es_field_operations.fieldSolve(es_field_state, integrated_charge);

  const mfem::GridFunction& potential = es_field_state.getPotential();
  const mfem::Vector expected_potential({0., -0.0414356, -0.0732684, -0.0957998, -0.109238, -0.113704,
    -0.109238, -0.0957998, -0.0732684, -0.0414356, 0.0});
  constexpr double relative_tolerance = 1.0e-3;
  ASSERT_NEAR(potential.DistanceTo(expected_potential) / expected_potential.Norml2(), 0.0, relative_tolerance);
}

// Mathematica
//
// solution :=
//  NDSolve[{-u''[x] + Exp[u[x]] == 2, u[0] == 0, u[1] == 0}, u[x], x]
// nodes := Subdivide[10]
// (u[x] /. solution)[[1]] /. x -> nodes
//
// {0., 0.0411263, 0.0726654, 0.0949505, 0.108224, 0.112632, 0.108224,
// 0.0949505, 0.0726654, 0.0411264, 6.75105*10^-9}
TEST(ElectrostaticFieldOperationsWithBoltzmannElectrons, ConstantChargeWithDirichletBoundaries) {
  constexpr int num_elems = 10;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);
  constexpr int hgrad_order = 1;
  Discretization es_discretization(&mesh, hgrad_order);
  std::unordered_map<std::string, int> side_name_to_boundary_attribute = getSideNameToBoundaryAttributeForInlineMeshes(
    mesh.Dimension());
  const int left_boundary_attribute = side_name_to_boundary_attribute["left"];
  const int right_boundary_attribute = side_name_to_boundary_attribute["right"];
  std::unordered_map<int, double> boundary_attribute_to_dirichlet_value{
    {left_boundary_attribute, 0.},
    {right_boundary_attribute, 0.}};
  auto dirichlet_bcs = std::make_unique<DirichletBoundaryConditionsConstant>(
    boundary_attribute_to_dirichlet_value, es_discretization);
  constexpr double reference_number_density = constants::permittivity / constants::elementary_charge;
  constexpr double temperature = constants::elementary_charge / constants::boltzmann_constant;
  ElectrostaticFieldOperationsWithBoltzmannElectrons es_field_operations(
    es_discretization,
    std::move(dirichlet_bcs),
    reference_number_density,
    temperature
  );
  ElectrostaticFieldState es_field_state(es_discretization);

  mfem::ConstantCoefficient charge_density(2. * constants::permittivity);
  mfem::LinearForm integrated_charge_linear_form(&es_discretization.getFeSpace());
  integrated_charge_linear_form.AddDomainIntegrator(new mfem::DomainLFIntegrator(charge_density));
  integrated_charge_linear_form.Assemble();
  IntegratedCharge integrated_charge(es_discretization);
  integrated_charge.setIntegratedCharge(integrated_charge_linear_form);

  es_field_operations.fieldSolve(es_field_state, integrated_charge);

  const mfem::GridFunction& potential = es_field_state.getPotential();
  const mfem::Vector expected_potential({0., 0.0411263, 0.0726654, 0.0949505, 0.108224, 0.112632, 0.108224,
    0.0949505, 0.0726654, 0.0411264, 0.0});
  constexpr double relative_tolerance = 1.0e-3;
  ASSERT_NEAR(potential.DistanceTo(expected_potential) / expected_potential.Norml2(), 0.0, relative_tolerance);
}
