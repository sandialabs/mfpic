#include <libmfpic/Constants.hpp>
#include <libmfpic/DGEulerAssembly.hpp>
#include <libmfpic/DGAssembly.hpp>
#include <libmfpic/DGEulerInitialConditionsFactory.hpp>
#include <libmfpic/DGEulerOperations.hpp>
#include <libmfpic/ElectrostaticFieldState.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/LowFidelityState.hpp>
#include <libmfpic/MeshFactory.hpp>
#include <libmfpic/SourcesFactory.hpp>
#include <libmfpic/Species.hpp>
#include <libmfpic/TestingUtils.hpp>

#include <gtest/gtest.h>
#include <mfem.hpp>

namespace {

using namespace mfpic;

TEST(DGEulerOperations, EulerChargeAssemblyWorksForConstantDensityOrder0In3D) {
  constexpr int dg_order = 0;
  constexpr int num_equations = 5;
  constexpr Species default_species{.charge = 2.0, .mass = 10.0, .specific_heat_ratio = 1.4};
  int nxyz = 3;

  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(nxyz, nxyz, nxyz, mfem::Element::HEXAHEDRON, 1, 1, 1);
  Discretization dg_discretization(&mesh, dg_order, FETypes::DG, num_equations);

  constexpr int charge_order = 1;
  Discretization charge_discretization(&mesh, charge_order, FETypes::HGRAD);

  mfem::FiniteElementSpace finite_element_space = dg_discretization.getFeSpace();
  std::shared_ptr<DGEulerAssembly> operator_ptr = std::make_shared<DGEulerAssembly>(finite_element_space, default_species);
  std::vector<std::shared_ptr<DGEulerAssembly>> dg_operators({operator_ptr});
  DGEulerOperations dg_euler_operations(charge_discretization, dg_operators);

  constexpr double number_density = 56.91;
  constexpr double temperature = 300;
  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters;
  list_of_parameters.push_back(std::make_unique<ConstantSourceParameters>(default_species, number_density, temperature));

  LowFidelityState current_state = buildEulerState(dg_discretization, list_of_parameters);

  IntegratedCharge charge_state = dg_euler_operations.assembleCharge(current_state);
  mfem::Vector integrated_charge_vector = charge_state.getIntegratedCharge();

  const int number_of_nodes = std::pow(nxyz + 1, 3);
  EXPECT_EQ(number_of_nodes, integrated_charge_vector.Size());

  const double expected_integrated_charge_interior = number_density * (1/std::pow(nxyz,3)) * default_species.charge;

  for (int i_node = 0; i_node < number_of_nodes; i_node++){
    const double* vertex_position = mesh.GetVertex(i_node);

    double expected_integrated_charge_at_node = expected_integrated_charge_interior;
    for (int i_dim = 0; i_dim < 3; ++i_dim) {
      if (vertex_position[i_dim] == 0. or vertex_position[i_dim] == 1.) {
        expected_integrated_charge_at_node *= 0.5;
      }
    }

    double integrated_charge_at_node = integrated_charge_vector[i_node];
    constexpr double relative_tolerance = 1e-15;
    EXPECT_NEAR_RELATIVE(expected_integrated_charge_at_node, integrated_charge_at_node, relative_tolerance);
  }
}

TEST(DGEulerOperations,EulerChargeAssemblyWorksForMultipleSpeciesConstantDensitiesOrder0In3D) {
  constexpr int dg_order = 0;
  constexpr int num_equations = 5;
  constexpr Species species_0{.charge = 2.0, .mass = 10.0};
  constexpr Species species_1{.charge = 532.0, .mass = 1920.0};
  constexpr int nxyz = 3;
  constexpr double length = 1.;
  constexpr double dx = length / nxyz;
  const double cell_volume = std::pow(dx, 3);
  const int num_nodes = std::pow(nxyz + 1, 3);

  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(nxyz, nxyz, nxyz, mfem::Element::HEXAHEDRON, length, length, length);
  Discretization dg_discretization(&mesh, dg_order, FETypes::DG, num_equations);

  constexpr int charge_order = 1;
  Discretization charge_discretization(&mesh, charge_order, FETypes::HGRAD);

  mfem::FiniteElementSpace finite_element_space = dg_discretization.getFeSpace();
  std::shared_ptr<DGEulerAssembly> operator_ptr_0 = std::make_shared<DGEulerAssembly>(finite_element_space, species_0);
  std::shared_ptr<DGEulerAssembly> operator_ptr_1 = std::make_shared<DGEulerAssembly>(finite_element_space, species_1);
  std::vector<std::shared_ptr<DGEulerAssembly>> dg_operators({operator_ptr_0, operator_ptr_1});
  DGEulerOperations dg_euler_operations(charge_discretization, dg_operators);

  constexpr double n_0 = 54.2342;
  constexpr double n_1 = 20.301;
  constexpr double temperature_0 = 320.1;
  constexpr double temperature_1 = 280.8;
  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters;
  list_of_parameters.push_back(std::make_unique<ConstantSourceParameters>(species_0, n_0, temperature_0));
  list_of_parameters.push_back(std::make_unique<ConstantSourceParameters>(species_1, n_1, temperature_1));
  LowFidelityState current_state = buildEulerState(dg_discretization, list_of_parameters);

  IntegratedCharge charge_state = dg_euler_operations.assembleCharge(current_state);
  mfem::Vector integrated_charge_vector = charge_state.getIntegratedCharge();
  EXPECT_EQ(num_nodes, integrated_charge_vector.Size());

  double expected_integrated_charge_interior = (n_0 * species_0.charge + n_1 * species_1.charge) * cell_volume;
  for (int i_node = 0; i_node < num_nodes; ++i_node) {
    const double* vertex_position = mesh.GetVertex(i_node);

    double expected_integrated_charge_at_node = expected_integrated_charge_interior;
    for (int i_dim = 0; i_dim < 3; ++i_dim) {
      if (vertex_position[i_dim] == 0. or vertex_position[i_dim] == length) {
        expected_integrated_charge_at_node *= 0.5;
      }
    }

    const double integrated_charge_at_node = integrated_charge_vector[i_node];
    constexpr double relative_tolerance = 1e-14;
    EXPECT_NEAR_RELATIVE(expected_integrated_charge_at_node, integrated_charge_at_node, relative_tolerance);
  }
}

TEST(DGEulerOperations,EulerChargeAssemblyWorksForLinearDensityOrder1In3D) {
  constexpr int dg_order = 1;
  constexpr int num_equations = 5;
  constexpr double specific_heat_ratio = 1.4;
  constexpr Species default_species{.charge = 1094.0, .mass = 9502.0, .specific_heat_ratio = specific_heat_ratio};

  constexpr int nxyz = 3;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(nxyz, nxyz, nxyz, mfem::Element::HEXAHEDRON, 1, 1, 1);
  Discretization dg_discretization(&mesh, dg_order, FETypes::DG, num_equations);

  constexpr int charge_order = 1;
  Discretization charge_discretization(&mesh, charge_order, FETypes::HGRAD);

  mfem::FiniteElementSpace finite_element_space = dg_discretization.getFeSpace();
  std::shared_ptr<DGEulerAssembly> operator_ptr = std::make_shared<DGEulerAssembly>(finite_element_space, default_species);
  std::vector<std::shared_ptr<DGEulerAssembly>> dg_operators({operator_ptr});
  DGEulerOperations dg_euler_operations(charge_discretization, dg_operators);

  constexpr double c0(53902.0), c1(923.0), c2(52.0), c3(892.52);
  auto solution_vec = [&](const mfem::Vector &x, mfem::Vector &y) {
    y = 0.;
    double rho = c0 + c1*x[0] + c2*x[1] + c3*x[2];
    y[0] = rho;
    y[4] = 10.43;
  };
  auto fluid_coeff = std::make_unique<mfem::VectorFunctionCoefficient>(num_equations, solution_vec);

  std::vector<std::pair<Species, std::unique_ptr<mfem::VectorCoefficient>>> list_of_species_and_coefficients;
  list_of_species_and_coefficients.push_back(std::make_pair(default_species, std::move(fluid_coeff)));
  LowFidelityState current_state(dg_discretization, list_of_species_and_coefficients);

  mfem::FiniteElementSpace test_finite_element_space = charge_discretization.getFeSpace();
  mfem::FunctionCoefficient test_rho([&](const mfem::Vector &x){
    return (default_species.charge / default_species.mass ) * (c0 + c1*x[0] + c2*x[1] + c3*x[2]);
  });
  mfem::LinearForm btest(&test_finite_element_space);
  btest.AddDomainIntegrator(new mfem::DomainLFIntegrator(test_rho));
  btest.Assemble();

  IntegratedCharge charge_state = dg_euler_operations.assembleCharge(current_state);
  mfem::Vector integrated_charge_vector = charge_state.getIntegratedCharge();
  integrated_charge_vector -= btest;
  constexpr double tolerance = 1e-12;
  EXPECT_NEAR(integrated_charge_vector.Norml2(), 0., tolerance);
}

TEST(DGEulerOperations,AccelerateDoesNotModifyInternalEnergy) {

  constexpr double tolerance = 1e-12;

  constexpr int dg_order = 1;
  constexpr int num_equations = 5;
  constexpr Species default_species{.charge = 2.0, .mass = 10.0};
  constexpr int nxyz = 3;

  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(nxyz, nxyz, nxyz, mfem::Element::HEXAHEDRON, 1, 1, 1);
  Discretization fluid_discretization(&mesh, dg_order, FETypes::DG, num_equations);

  constexpr int potential_order = 1;
  Discretization potential_discretization(&mesh, potential_order, FETypes::HGRAD);

  auto operator_ptr = std::make_shared<DGEulerAssembly>(fluid_discretization.getFeSpace(), default_species);
  std::vector<std::shared_ptr<DGEulerAssembly>> dg_operators({operator_ptr});
  DGEulerOperations dg_euler_operations(potential_discretization, dg_operators);

  constexpr mfem::real_t dx = -3.4;
  constexpr mfem::real_t dy = 8.2;
  constexpr mfem::real_t dz = 0.8;
  auto linear_func = [&](const mfem::Vector &x){ return dx * x[0] + dy * x[1] + dz * x[2]; };
  mfem::FunctionCoefficient potential_coeff(linear_func);
  mfem::GridFunction potential(&potential_discretization.getFeSpace());
  potential.ProjectCoefficient(potential_coeff);
  mfem::Vector e{-dx, -dy, -dz};

  ElectrostaticFieldState es_field_state(potential_discretization);
  es_field_state.setPotential(potential);

  auto compute_internal_energy_at_quad_pts = [&](const mfem::GridFunction & fluid_grid_function) {
    std::vector<double> internal_energy;

    mfem::Array<int> vector_dofs;
    mfem::DenseMatrix ip_locations, fluid_vals_at_ips;
    const auto & fe_space = fluid_discretization.getFeSpace();
    for (int element=0; element < fe_space.GetNE(); ++element) {
      fe_space.GetElementVDofs(element, vector_dofs);
      const mfem::IntegrationRule &integration_rule = mfem::IntRules.Get(fe_space.GetFE(element)->GetGeomType(), 2*fe_space.GetFE(element)->GetOrder());
      mfem::ElementTransformation* element_transformation = fe_space.GetElementTransformation(element);
      element_transformation->Transform(integration_rule, ip_locations);
      fluid_grid_function.GetVectorValues(*element_transformation, integration_rule, fluid_vals_at_ips);

      mfem::Vector fluid_vals_at_single_pt(5);

      for (int ipoint = 0; ipoint < integration_rule.GetNPoints(); ++ipoint) {
        fluid_vals_at_ips.GetColumn(ipoint, fluid_vals_at_single_pt);
        const double density = fluid_vals_at_single_pt[0];
        const mfem::Vector momentum {fluid_vals_at_single_pt[1], fluid_vals_at_single_pt[2], fluid_vals_at_single_pt[3]};
        const mfem::Vector velocity {momentum[0]/density, momentum[1]/density, momentum[2]/density};
        const double ke = .5 * (momentum * velocity);
        const double ie = fluid_vals_at_single_pt[4] - ke;
        internal_energy.push_back(ie);
      }
    }
    return internal_energy;
  };

  constexpr mfem::real_t c0(12.7), c1(-9.4), c2(2.2), c3(9.1);
  auto linear_vec = [&](const mfem::Vector &x, mfem::Vector &y) {
    mfem::real_t base_val = c0 + c1 * x[0] + c2 * x[1] + c3 * x[2];
    for (int i = 0; i < num_equations; ++i) {
      y[i] = (i+1) * base_val;
    }
  };
  auto fluid_coeff = std::make_unique<mfem::VectorFunctionCoefficient>(num_equations,linear_vec);

  std::vector<std::pair<Species, std::unique_ptr<mfem::VectorCoefficient>>> list_of_species_and_coefficients;
  list_of_species_and_coefficients.push_back(std::make_pair(default_species, std::move(fluid_coeff)));
  LowFidelityState state_before(fluid_discretization, list_of_species_and_coefficients);
  LowFidelitySpeciesState& species_state_before = state_before.getSpeciesState(0);
  mfem::GridFunction fluid_grid_function = species_state_before.getGridFunction();

  const std::vector<double> internal_energy_before = compute_internal_energy_at_quad_pts(fluid_grid_function);

  LowFidelityState state_after = dg_euler_operations.accelerate(1.0, state_before, es_field_state);
  const LowFidelitySpeciesState& species_state_after = state_after.getSpeciesState(0);
  const mfem::GridFunction& final_grid_function = species_state_after.getGridFunction();

  mfem::Vector diff = final_grid_function;
  diff -= fluid_grid_function;
  const double diff_l2 = diff.Norml2();

  // the coefficients have changed
  EXPECT_GE(diff_l2, 1e-1);

  const auto internal_energy_after = compute_internal_energy_at_quad_pts(final_grid_function);

  EXPECT_EQ(internal_energy_before.size(), internal_energy_after.size());

  for (size_t i = 0; i < internal_energy_before.size(); ++i) {
    EXPECT_NEAR(internal_energy_before[i], internal_energy_after[i], tolerance);
  }

}

TEST(DGEulerOperations, MoveConstant) {
  MeshParameters mesh_parameters{
    .mesh_type = "hex",
    .lengths = {1., 1., 1.},
    .num_elements = {3, 3, 3},
    .periodic_dims = {0, 1, 2}};

  mfem::Mesh mesh = buildMesh(mesh_parameters);

  constexpr int dg_order = 0;
  const int num_equations = euler::ConservativeVariables::NUM_VARS;

  Discretization dg_discretization(&mesh, dg_order, FETypes::DG, num_equations);

  Species electron_species{.charge = -constants::elementary_charge, .mass = constants::electron_mass};
  constexpr double number_density = 1e22;
  constexpr double temperature = 300;
  const mfem::Vector bulk_velocity{1., 2., 3.};
  std::vector<std::unique_ptr<SourceParameters>> list_of_ic_parameters;
  list_of_ic_parameters.push_back(
    std::make_unique<ConstantSourceParameters>(electron_species, number_density, temperature, bulk_velocity));
  const LowFidelityState dg_euler_state = buildEulerState(dg_discretization, list_of_ic_parameters);
  const LowFidelitySpeciesState& species_state = dg_euler_state.getSpeciesState(0);
  const mfem::GridFunction& grid_function = species_state.getGridFunction();

  const int charge_order = 1;
  Discretization charge_discretization(&mesh, charge_order, FETypes::HGRAD);

  auto dg_euler_assembly = std::make_shared<DGEulerAssembly>(dg_discretization.getFeSpace(), electron_species);
  std::vector<std::shared_ptr<DGEulerAssembly>> dg_assemblers({dg_euler_assembly});
  DGEulerOperations dg_euler_operations(charge_discretization, dg_assemblers);

  constexpr double dt = 1e-5;
  const LowFidelityState dg_euler_state_moved = dg_euler_operations.move(dt, dg_euler_state);
  const LowFidelitySpeciesState& species_state_moved = dg_euler_state_moved.getSpeciesState(0);
  const mfem::GridFunction& grid_function_moved = species_state_moved.getGridFunction();

  mfem::VectorGridFunctionCoefficient expected_solution(&grid_function);
  const double error = grid_function_moved.ComputeL2Error(expected_solution);
  constexpr double tolerance = 1e-16;
  EXPECT_LE(error, tolerance);
}

} // namespace
