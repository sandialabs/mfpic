#include <libmfpic/Constants.hpp>
#include <libmfpic/DirichletBoundaryConditionsConstant.hpp>
#include <libmfpic/ElectrostaticFieldOperations.hpp>
#include <libmfpic/ElectrostaticFieldState.hpp>
#include <libmfpic/IntegratedCharge.hpp>
#include <libmfpic/MeshFactory.hpp>
#include <libmfpic/MeshUtilities.hpp>
#include <libmfpic/Pinning.hpp>
#include <libmfpic/RandomNumberGenerator.hpp>

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

mfem::Vector excludeBoundaryEntries(const mfem::Vector& full_vector, mfem::Array<int>& boundary_dofs) {
  mfem::Array<int> marker(full_vector.Size());
  marker = 0;
  for (int i = 0; i < boundary_dofs.Size(); i++)
  {
    marker[boundary_dofs[i]] = 1;
  }

  mfem::Array<int> interior_dofs;
  for (int i = 0; i < full_vector.Size(); i++)
  {
      if (!marker[i])
      {
        interior_dofs.Append(i);
      }
  }

  mfem::Vector interior_values(interior_dofs.Size());

  for (int i = 0; i < interior_dofs.Size(); i++)
    interior_values(i) = full_vector(interior_dofs[i]);

  return interior_values;
}

// 2 eps charge
// zero on boundaries
// potential should be -(x - 1/2)^2 + 1/4
TEST(ElectrostaticFieldOperations, NoChargeErrorAfterPoissonSolveDirichlet) {
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
  auto dirichlet_dofs = dirichlet_bcs->getDirichletBoundaryDofIndices();

  ElectrostaticFieldOperations es_field_operations(es_discretization, std::move(dirichlet_bcs));
  ElectrostaticFieldState es_field_state(es_discretization);

  mfem::ConstantCoefficient charge_density(2. * constants::permittivity);
  mfem::LinearForm integrated_charge_linear_form(&es_discretization.getFeSpace());
  integrated_charge_linear_form.AddDomainIntegrator(new mfem::DomainLFIntegrator(charge_density));
  integrated_charge_linear_form.Assemble();

  IntegratedCharge integrated_charge(es_discretization);
  integrated_charge.setIntegratedCharge(integrated_charge_linear_form);

  es_field_operations.fieldSolve(es_field_state, integrated_charge);

  mfem::Vector integrated_ghost_charge = es_field_operations.computeIntegratedGhostCharge(es_field_state, integrated_charge);
  const mfem::Vector integrated_ghost_charge_interior = excludeBoundaryEntries(integrated_ghost_charge, dirichlet_dofs);

  const double charge_error_l2_norm = integrated_ghost_charge_interior.Norml2();
  const double integrated_charge_l2_norm = integrated_charge_linear_form.Norml2();
  const double relative_error = charge_error_l2_norm / integrated_charge_l2_norm;

  constexpr double relative_tolerance = 1e-11;
  EXPECT_LT(abs(relative_error), relative_tolerance);
}

// 2 eps charge
// zero on boundaries
// potential should be -(x - 1/2)^2 + 1/4
// we then modify the potential by 2x
TEST(ElectrostaticFieldOperations, QuadraticOffsetGivesConstantChargeError) {
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
  auto dirichlet_dofs = dirichlet_bcs->getDirichletBoundaryDofIndices();

  ElectrostaticFieldOperations es_field_operations(es_discretization, std::move(dirichlet_bcs));
  ElectrostaticFieldState es_field_state(es_discretization);

  mfem::ConstantCoefficient charge_density(2. * constants::permittivity);
  mfem::LinearForm integrated_charge_linear_form(&es_discretization.getFeSpace());
  integrated_charge_linear_form.AddDomainIntegrator(new mfem::DomainLFIntegrator(charge_density));
  integrated_charge_linear_form.Assemble();

  IntegratedCharge integrated_charge(es_discretization);
  integrated_charge.setIntegratedCharge(integrated_charge_linear_form);

  es_field_operations.fieldSolve(es_field_state, integrated_charge);

  mfem::GridFunction& potential = es_field_state.getPotential();
  potential *= 2.;

  mfem::Vector integrated_ghost_charge = es_field_operations.computeIntegratedGhostCharge(es_field_state, integrated_charge);

  // we doubled the potential so should have an error equal to the charge we started with
  const mfem::Vector expected_ghost_charge = integrated_charge.getIntegratedCharge();
  const mfem::Vector integrated_ghost_charge_interior = excludeBoundaryEntries(integrated_ghost_charge, dirichlet_dofs);
  const mfem::Vector expected_ghost_charge_interior = excludeBoundaryEntries(expected_ghost_charge, dirichlet_dofs);

  constexpr double solver_tolerance = 1e-12;
  constexpr double pointwise_tolerance = 10 * solver_tolerance;
  for (int i = 0; i < integrated_ghost_charge_interior.Size(); ++i) {
    EXPECT_LT(integrated_ghost_charge_interior(i) - expected_ghost_charge_interior(i), pointwise_tolerance * abs(expected_ghost_charge_interior(i)));
  }

}

TEST(ElectrostaticFieldOperations, NoChargeErrorAfterPoissonSolvePeriodic) {
  MeshParameters mesh_parameters{
    .mesh_type = "line",
    .lengths = {1.},
    .num_elements = {15},
    .periodic_dims = {0}};

  mfem::Mesh mesh = buildMesh(mesh_parameters);

  constexpr int hgrad_order = 1;
  Discretization es_discretization(&mesh, hgrad_order);

  auto pinning = std::make_unique<Pinning>();
  ElectrostaticFieldOperations es_field_operations(es_discretization, std::move(pinning));
  ElectrostaticFieldState es_field_state(es_discretization);

  mfem::FunctionCoefficient charge_density([](const mfem::Vector& x){
    return 4. * M_PI * M_PI * constants::permittivity * cos(2. * M_PI * x[0]);
  });
  mfem::LinearForm integrated_charge_linear_form(&es_discretization.getFeSpace());
  integrated_charge_linear_form.AddDomainIntegrator(new mfem::DomainLFIntegrator(charge_density));
  integrated_charge_linear_form.Assemble();

  IntegratedCharge integrated_charge(es_discretization);
  integrated_charge.setIntegratedCharge(integrated_charge_linear_form);

  es_field_operations.fieldSolve(es_field_state, integrated_charge);

  mfem::Vector integrated_ghost_charge = es_field_operations.computeIntegratedGhostCharge(es_field_state, integrated_charge);

  mfem::Vector integrated_charge_vector = integrated_charge.getIntegratedCharge();
  constexpr double solver_tolerance = 1e-12;
  constexpr double pointwise_tolerance = 10 * solver_tolerance;
  for (int i = 0; i < integrated_ghost_charge.Size(); ++i) {
    EXPECT_LT(abs(integrated_ghost_charge(i)), pointwise_tolerance * abs(integrated_charge_vector(i)));
  }

  mfem::GridFunction ghost_charge_density = es_field_operations.computeGhostChargeDensity(es_field_state, integrated_charge);
  const double ghost_charge_density_norm = ghost_charge_density.Norml2();
  constexpr double tolerance = 1e-16;
  EXPECT_LT(ghost_charge_density_norm, tolerance);
}

TEST(ElectrostaticFieldOperations, FieldSolveCorrectForPeriodicWithPinning) {
  MeshParameters mesh_parameters{
    .mesh_type = "line",
    .lengths = {1.},
    .num_elements = {60},
    .periodic_dims = {0}};

  mfem::Mesh mesh = buildMesh(mesh_parameters);

  constexpr int hgrad_order = 1;
  Discretization es_discretization(&mesh, hgrad_order);

  auto pinning = std::make_unique<Pinning>();
  ElectrostaticFieldOperations es_field_operations(es_discretization, std::move(pinning));
  ElectrostaticFieldState es_field_state(es_discretization);

  mfem::FunctionCoefficient charge_density([](const mfem::Vector& x){
    return 4. * M_PI * M_PI * constants::permittivity * cos(2. * M_PI * x[0]);
  });
  mfem::LinearForm integrated_charge_linear_form(&es_discretization.getFeSpace());
  integrated_charge_linear_form.AddDomainIntegrator(new mfem::DomainLFIntegrator(charge_density));
  integrated_charge_linear_form.Assemble();

  IntegratedCharge integrated_charge(es_discretization);
  integrated_charge.setIntegratedCharge(integrated_charge_linear_form);

  es_field_operations.fieldSolve(es_field_state, integrated_charge);

  mfem::FunctionCoefficient exact_potential([](const mfem::Vector& x){
    const double x_pin = 0.;
    const double pin_value = cos(2. * M_PI * x_pin);
    return cos(2. * M_PI * x[0]) - pin_value;
  });

  const mfem::GridFunction& potential = es_field_state.getPotential();
  const double l2_error = potential.ComputeL2Error(exact_potential);
  constexpr double expected_error = 0.;
  constexpr double absolute_tolerance = 1e-3;
  EXPECT_NEAR(l2_error, expected_error, absolute_tolerance);

}

TEST(ElectrostaticFieldOperations, CheckCompatibilityEnforcementIsAProjection) {

  MeshParameters mesh_parameters{
    .mesh_type = "line",
    .lengths = {1.},
    .num_elements = {60},
    .periodic_dims = {0}};

  mfem::Mesh mesh = buildMesh(mesh_parameters);

  constexpr int hgrad_order = 1;
  Discretization es_discretization(&mesh, hgrad_order);
  auto pinning = std::make_unique<Pinning>();
  ElectrostaticFieldOperations es_field_operations(es_discretization, std::move(pinning));

  IntegratedCharge integrated_charge(es_discretization);
  mfem::Vector integrated_charge_vector = integrated_charge.getIntegratedCharge();
  std::random_device rd;
  RandomNumberGenerator gen(rd());
  std::uniform_real_distribution<> dis(1.0, 2.0);
  for (int i = 0; i < integrated_charge_vector.Size(); ++i) {
    integrated_charge_vector[i] = dis(gen);
  }

  // check this is a projection, i.e., P^2 = P
  es_field_operations.enforceCompatibilityOnIntegratedCharge(integrated_charge_vector);
  const auto first_time = integrated_charge_vector;
  es_field_operations.enforceCompatibilityOnIntegratedCharge(integrated_charge_vector);
  mfem::Vector diff = first_time;
  diff.Add(-1,integrated_charge_vector);
  EXPECT_NEAR(diff.Norml2(),0.,1e-14);

}

TEST(ElectrostaticFieldOperations, CheckNullSpaceVectorIsOne) {

  MeshParameters mesh_parameters{
    .mesh_type = "line",
    .lengths = {1.},
    .num_elements = {60},
    .periodic_dims = {0}};

  mfem::Mesh mesh = buildMesh(mesh_parameters);

  constexpr int hgrad_order = 1;
  Discretization es_discretization(&mesh, hgrad_order);
  auto pinning = std::make_unique<Pinning>();
  ElectrostaticFieldOperations es_field_operations(es_discretization, std::move(pinning));

  mfem::GridFunction potential(&es_discretization.getFeSpace());
  mfem::ConstantCoefficient one(1.0);
  potential.ProjectCoefficient(one);

  ElectrostaticFieldState es_field_state(es_discretization);
  es_field_state.setPotential(potential);

  // this is the energy norm so it being zero implies we're in the null space
  const double field_energy = es_field_operations.fieldEnergy(es_field_state);
  EXPECT_NEAR(field_energy, 0., 1e-14);

  // after correction should be zero if we started in the null space
  es_field_operations.enforceCompatibilityOnIntegratedCharge(potential);
  EXPECT_NEAR(potential.Norml2(), 0., 1e-14);

}

}
