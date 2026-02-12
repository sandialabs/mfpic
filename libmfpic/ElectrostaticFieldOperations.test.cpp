#include <libmfpic/Constants.hpp>
#include <libmfpic/DirichletBoundaryConditionsConstant.hpp>
#include <libmfpic/ElectrostaticFieldOperations.hpp>
#include <libmfpic/ElectrostaticFieldState.hpp>
#include <libmfpic/MeshUtilities.hpp>
#include <libmfpic/Pinning.hpp>

#include <mfem/mfem.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

TEST(ElectrostaticFieldOperations, ZeroChargeWithNaturalBoundariesGivesZeroPotential) {
  constexpr int num_elems = 10;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);

  constexpr int hgrad_order = 1;
  Discretization es_discretization(&mesh, hgrad_order);

  auto pinning = std::make_unique<Pinning>();
  ElectrostaticFieldOperations es_field_operations(es_discretization, std::move(pinning));
  ElectrostaticFieldState es_field_state(es_discretization);

  IntegratedCharge integrated_charge(es_discretization);
  integrated_charge.setIntegratedChargeValue(0.);

  es_field_operations.fieldSolve(es_field_state, integrated_charge);

  mfem::GridFunction potential = es_field_state.getPotential();
  mfem::ConstantCoefficient exact_potential(0.);
  const double l2_error = potential.ComputeL2Error(exact_potential);

  constexpr double expected_error = 0.;
  EXPECT_EQ(expected_error, l2_error);
}

// if charge density = 0 and potential at x = 0 is l and potential at x = 1 is r
// then exact solution is potential = l (1 - x) + r x
TEST(ElectrostaticFieldOperations, ZeroChargeWithDirichletBoundariesGivesLinearPotential) {
  constexpr int num_elems = 20;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);

  constexpr int hgrad_order = 1;
  Discretization es_discretization(&mesh, hgrad_order);

  std::unordered_map<std::string, int> side_name_to_boundary_attribute = getSideNameToBoundaryAttributeForInlineMeshes(
    mesh.Dimension());

  const int left_boundary_attribute = side_name_to_boundary_attribute["left"];
  const int right_boundary_attribute = side_name_to_boundary_attribute["right"];
  const double left_boundary_value = 2.1;
  const double right_boundary_value = 3.4;
  std::unordered_map<int, double> boundary_attribute_to_dirichlet_value{
    {left_boundary_attribute, left_boundary_value},
    {right_boundary_attribute, right_boundary_value}};
  auto dirichlet_bcs = std::make_unique<DirichletBoundaryConditionsConstant>(
    boundary_attribute_to_dirichlet_value, es_discretization);

  ElectrostaticFieldOperations es_field_operations(es_discretization, std::move(dirichlet_bcs));
  ElectrostaticFieldState es_field_state(es_discretization);

  IntegratedCharge integrated_charge(es_discretization);
  integrated_charge.setIntegratedChargeValue(0.);

  es_field_operations.fieldSolve(es_field_state, integrated_charge);

  mfem::GridFunction potential = es_field_state.getPotential();

  mfem::FunctionCoefficient exact_potential([left_boundary_value, right_boundary_value](const mfem::Vector& x) -> mfem::real_t {
    return left_boundary_value * (1 - x[0]) + right_boundary_value * x[0];
  });
  const double l2_error = potential.ComputeL2Error(exact_potential);

  constexpr double tolerance = 1e-12;
  EXPECT_LT(l2_error, tolerance);
}

// 2 eps charge
// zero on boundaries
// potential should be -(x - 1/2)^2 + 1/4
TEST(ElectrostaticFieldOperations, ConstantChargeGivesQuadraticPotential) {
  constexpr int num_elems = 20;
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

  ElectrostaticFieldOperations es_field_operations(es_discretization, std::move(dirichlet_bcs));
  ElectrostaticFieldState es_field_state(es_discretization);

  mfem::ConstantCoefficient charge_density(2. * constants::permittivity);
  mfem::LinearForm integrated_charge_linear_form(&es_discretization.getFeSpace());
  integrated_charge_linear_form.AddDomainIntegrator(new mfem::DomainLFIntegrator(charge_density));
  integrated_charge_linear_form.Assemble();

  IntegratedCharge integrated_charge(es_discretization);
  integrated_charge.setIntegratedCharge(integrated_charge_linear_form);

  es_field_operations.fieldSolve(es_field_state, integrated_charge);

  mfem::GridFunction potential = es_field_state.getPotential();
  mfem::FunctionCoefficient exact_potential([](const mfem::Vector& x) -> mfem::real_t {
    return -1. * (x[0] - 0.5) * (x[0] - 0.5) + 0.25;
  });
  const double l2_error = potential.ComputeL2Error(exact_potential);

  constexpr double tolerance = 0.001;
  EXPECT_LT(l2_error, tolerance);

  mfem::Vector potential_nodal_values;
  potential.GetNodalValues(potential_nodal_values);

  mfem::GridFunction exact_potential_grid_function(&es_discretization.getFeSpace());
  exact_potential_grid_function.ProjectCoefficient(exact_potential);

  mfem::Vector exact_potential_nodal_values;
  exact_potential_grid_function.GetNodalValues(exact_potential_nodal_values);

  mfem::Vector error_nodal_values = potential_nodal_values;
  error_nodal_values -= exact_potential_nodal_values;

  const double nodal_error_l2_norm = error_nodal_values.Norml2();
  constexpr double expected_nodal_error_l2_norm = 0.;
  constexpr double absolute_tolerance = 1e-13;
  EXPECT_NEAR(expected_nodal_error_l2_norm, nodal_error_l2_norm, absolute_tolerance);
}

TEST(ElectrostaticFieldOperations, ZeroPotentialHasZeroFieldEnergy) {
  constexpr int num_elems = 4;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);
  constexpr int hgrad_order = 1;
  Discretization es_discretization(&mesh, hgrad_order);
  std::unique_ptr<DirichletBoundaryConditionsConstant> empty_dirichlet_bcs;
  ElectrostaticFieldOperations es_field_operations(es_discretization, std::move(empty_dirichlet_bcs));
  ElectrostaticFieldState es_field_state(es_discretization);

  es_field_state.getPotential() = 0.0;
  const double field_energy = es_field_operations.fieldEnergy(es_field_state);

  ASSERT_DOUBLE_EQ(field_energy, 0.0);
}

TEST(ElectrostaticFieldOperations, ConstantPotentialHasZeroFieldEnergy) {
  constexpr int num_elems = 4;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);
  constexpr int hgrad_order = 1;
  Discretization es_discretization(&mesh, hgrad_order);
  std::unique_ptr<DirichletBoundaryConditionsConstant> empty_dirichlet_bcs;
  ElectrostaticFieldOperations es_field_operations(es_discretization, std::move(empty_dirichlet_bcs));
  ElectrostaticFieldState es_field_state(es_discretization);

  es_field_state.getPotential() = 505.0;
  const double field_energy = es_field_operations.fieldEnergy(es_field_state);

  constexpr double absolute_tolerance = 1e-15;
  ASSERT_NEAR(field_energy, 0.0, absolute_tolerance);
}

TEST(ElectrostaticFieldOperations, LinearPotentialFieldEnergy) {
  constexpr int num_elems = 4;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);
  constexpr int hgrad_order = 1;
  Discretization es_discretization(&mesh, hgrad_order);
  std::unique_ptr<DirichletBoundaryConditionsConstant> empty_dirichlet_bcs;
  ElectrostaticFieldOperations es_field_operations(es_discretization, std::move(empty_dirichlet_bcs));
  ElectrostaticFieldState es_field_state(es_discretization);

  mfem::GridFunction& potential = es_field_state.getPotential();
  constexpr double gradient = 1000.0;
  mfem::FunctionCoefficient linear_function([=](const mfem::Vector& position) -> double { return gradient * position[0]; });
  potential.ProjectCoefficient(linear_function);
  const double field_energy = es_field_operations.fieldEnergy(es_field_state);

  constexpr double domain_length = 1.0;
  constexpr double expected_field_energy = 0.5 * constants::permittivity * gradient * gradient * domain_length;
  ASSERT_DOUBLE_EQ(field_energy, expected_field_energy);
}

}
