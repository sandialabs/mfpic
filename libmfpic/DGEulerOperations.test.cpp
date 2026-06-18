#include <libmfpic/Constants.hpp>
#include <libmfpic/DGEulerAssembly.hpp>
#include <libmfpic/DGAssembly.hpp>
#include <libmfpic/DGEulerInitialConditionsFactory.hpp>
#include <libmfpic/DGEulerOperations.hpp>
#include <libmfpic/DGEulerOperationsFactory.hpp>
#include <libmfpic/ElectrostaticFieldState.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/LowFidelityState.hpp>
#include <libmfpic/MeshFactory.hpp>
#include <libmfpic/SourcesFactory.hpp>
#include <libmfpic/Species.hpp>
#include <libmfpic/TestingUtils.hpp>

#include <gtest/gtest.h>
#include <mfem.hpp>

#include <ranges>

namespace {

using namespace mfpic;

Species electron_species{.charge = -constants::elementary_charge, .mass = constants::electron_mass};
Species proton_species{.charge = constants::elementary_charge, .mass = constants::proton_mass};

std::vector<std::pair<Species, std::unique_ptr<mfem::VectorCoefficient>>> empty_sources{};

TEST(DGEulerOperations, EulerChargeAssemblyWorksForConstantDensityOrder0In3D) {
  constexpr int dg_order = 0;
  constexpr int num_equations = 5;
  const Species default_species{.charge = 2.0, .mass = 10.0, .specific_heat_ratio = 1.4};
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

  LowFidelityState state = buildEulerState(dg_discretization, list_of_parameters);

  IntegratedCharge charge_state = dg_euler_operations.assembleCharge(state);
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
  const Species species_0{.charge = 2.0, .mass = 10.0};
  const Species species_1{.charge = 532.0, .mass = 1920.0};
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
  LowFidelityState state = buildEulerState(dg_discretization, list_of_parameters);

  IntegratedCharge charge_state = dg_euler_operations.assembleCharge(state);
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
  const Species default_species{.charge = 1094.0, .mass = 9502.0, .specific_heat_ratio = specific_heat_ratio};

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
  LowFidelityState state(dg_discretization, list_of_species_and_coefficients);

  mfem::FiniteElementSpace test_finite_element_space = charge_discretization.getFeSpace();
  mfem::FunctionCoefficient test_rho([&](const mfem::Vector &x){
    return (default_species.charge / default_species.mass ) * (c0 + c1*x[0] + c2*x[1] + c3*x[2]);
  });
  mfem::LinearForm btest(&test_finite_element_space);
  btest.AddDomainIntegrator(new mfem::DomainLFIntegrator(test_rho));
  btest.Assemble();

  IntegratedCharge charge_state = dg_euler_operations.assembleCharge(state);
  mfem::Vector integrated_charge_vector = charge_state.getIntegratedCharge();
  integrated_charge_vector -= btest;
  constexpr double tolerance = 1e-12;
  EXPECT_NEAR(integrated_charge_vector.Norml2(), 0., tolerance);
}

TEST(DGEulerOperations,AccelerateDoesNotModifyInternalEnergy) {

  constexpr double tolerance = 1e-12;

  constexpr int dg_order = 1;
  constexpr int num_equations = 5;
  const Species default_species{.charge = 2.0, .mass = 10.0};
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

TEST(DGEulerOperations, MoveAccelerate_ZeroFieldsReproducesMove) {
  constexpr double length = 1.5;
  MeshParameters mesh_parameters{.mesh_type = "line", .lengths = {length}, .num_elements = {20}, .periodic_dims = {0}};
  mfem::Mesh mesh = buildMesh(mesh_parameters);

  constexpr int dg_order = 0;
  const int num_equations = euler::ConservativeVariables::NUM_VARS;
  Discretization dg_discretization(&mesh, dg_order, FETypes::DG, num_equations);

  const mfem::Vector center{0.5 * length};
  constexpr double standard_deviation = 0.1 * length;
  const SourceStateParameters offsets{.number_density = 1.8e21, .bulk_velocity = mfem::Vector({10.2, 0., 0.})};
  const SourceStateParameters heights{.number_density = 9.2e19, .bulk_velocity = mfem::Vector({4.8, 0., 0.})};
  constexpr double pressure_offset = 29713.;
  constexpr double pressure_height = 3479.;
  std::vector<std::unique_ptr<SourceParameters>> list_of_ic_parameters;
  list_of_ic_parameters.push_back(
    std::make_unique<GaussianSourceParameters>(
      electron_species, center, standard_deviation, offsets, heights, pressure_offset, pressure_height));

  const LowFidelityState initial_state = buildEulerState(dg_discretization, list_of_ic_parameters);

  constexpr int es_order = 1;
  Discretization es_discretization(&mesh, es_order, FETypes::HGRAD);
  const std::vector<Species> species_list = initial_state.getSpeciesList();
  std::vector<std::vector<std::unique_ptr<DGBC>>> empty_bcs(species_list.size());
  std::unique_ptr<LowFidelityOperations> dg_euler_operations = buildDGEulerOperations(
    dg_discretization, es_discretization, species_list, empty_bcs, empty_sources);

  ElectrostaticFieldState es_field_state_zero(es_discretization);
  constexpr double dt = 1e-8;

  const LowFidelityState state_moved = dg_euler_operations->move(dt, initial_state);
  const LowFidelityState state_moved_accelerated = dg_euler_operations->moveAccelerate(dt, initial_state, es_field_state_zero);

  const LowFidelitySpeciesState& species_state_moved = state_moved.getSpeciesState(0);
  const LowFidelitySpeciesState& species_state_moved_accelerated = state_moved_accelerated.getSpeciesState(0);
  const mfem::GridFunction& grid_function_moved = species_state_moved.getGridFunction();
  const mfem::GridFunction& grid_function_moved_acclerated = species_state_moved_accelerated.getGridFunction();
  for (int i = 0; i < grid_function_moved.Size(); ++i) {
    EXPECT_DOUBLE_EQ(grid_function_moved[i], grid_function_moved_acclerated[i]);
  }
}

class AppliedElectricField : public ElectromagneticFieldsEvaluator {
public:
  AppliedElectricField(const mfem::Vector e_field) : e_field_(e_field) {}

  virtual mfem::Vector getEFieldAt(const mfem::Vector& /*position*/, const int /*element_index*/) const override {
    return e_field_;
  }

  virtual mfem::Vector getBFieldAt(const mfem::Vector& /*position*/, const int /*element_index*/) const override {
    mfem::Vector b_field(3);
    b_field = 0.;
    return b_field;
  }
private:
  mfem::Vector e_field_;
};

TEST(DGEulerOperations, MoveAccelerate_ConstantEFieldAcceleratesConstantFluidCorrectly) {
  constexpr double length = 1.5;
  MeshParameters mesh_parameters{.mesh_type = "line", .lengths = {length}, .num_elements = {20}, .periodic_dims = {0}};
  mfem::Mesh mesh = buildMesh(mesh_parameters);

  constexpr int dg_order = 0;
  const int num_equations = euler::ConservativeVariables::NUM_VARS;
  Discretization dg_discretization(&mesh, dg_order, FETypes::DG, num_equations);

  constexpr double number_density = 3.4e23;
  constexpr double temperature = 321;
  const mfem::Vector bulk_velocity{8.4e2, 1.9e3, -5.7e4};
  std::vector<std::unique_ptr<SourceParameters>> list_of_ic_parameters;
  list_of_ic_parameters.push_back(
    std::make_unique<ConstantSourceParameters>(electron_species, number_density, temperature, bulk_velocity));

  const LowFidelityState initial_state = buildEulerState(dg_discretization, list_of_ic_parameters);

  constexpr int es_order = 1;
  Discretization es_discretization(&mesh, es_order, FETypes::HGRAD);
  const std::vector<Species> species_list = initial_state.getSpeciesList();
  std::vector<std::vector<std::unique_ptr<DGBC>>> empty_bcs(species_list.size());
  std::unique_ptr<LowFidelityOperations> dg_euler_operations = buildDGEulerOperations(
    dg_discretization, es_discretization, species_list, empty_bcs, empty_sources);

  const double plasma_frequency = sqrt(
    (number_density * electron_species.charge * electron_species.charge) / (electron_species.mass * constants::permittivity));
  const double plasma_period = 2.0 * M_PI / plasma_frequency;
  const double dt = 0.1 * plasma_period;

  const mfem::Vector e_field{2.3e3, 4.5e4, 6.7e4};
  const AppliedElectricField applied_e_field(e_field);
  const LowFidelityState state_moved_accelerated = dg_euler_operations->moveAccelerate(dt, initial_state, applied_e_field);

  const LowFidelitySpeciesState& species_state_moved_accelerated = state_moved_accelerated.getSpeciesState(0);
  const mfem::GridFunction& grid_function_moved_acclerated = species_state_moved_accelerated.getGridFunction();

  mfem::Vector primitive_state = euler::constructPrimitiveState(number_density, bulk_velocity, temperature);
  mfem::Vector conservative_state = euler::convertFromPrimitiveToConservative(primitive_state, electron_species);
  const mfem::Vector initial_momentum_density = euler::getMomentumDensityFromConservativeState(conservative_state);
  const double initial_total_energy_density = conservative_state[euler::ConservativeVariables::TOTAL_ENERGY_DENSITY];

  const double expected_mass_density = conservative_state[euler::ConservativeVariables::MASS_DENSITY];
  mfem::Vector expected_momentum_density = initial_momentum_density;
  expected_momentum_density.Add(dt * electron_species.charge * number_density, e_field);
  const double charge_to_mass = electron_species.charge / electron_species.mass;
  const double total_energy_density_increment = dt * charge_to_mass * (initial_momentum_density * e_field);
  const double expected_total_energy_density = initial_total_energy_density + total_energy_density_increment;

  const mfem::Vector expected_conservative_state = euler::constructConservativeState(
    expected_mass_density, expected_momentum_density, expected_total_energy_density);
  mfem::VectorConstantCoefficient expected_coefficient(expected_conservative_state);
  const double error = grid_function_moved_acclerated.ComputeL2Error(expected_coefficient);
  constexpr double tolerance = 1e-16;
  EXPECT_LT(error, tolerance);
}

TEST(DGEulerOperations, computeRHS_ConstantStateWithNoFieldsGivesZeroRHS) {
  constexpr double length = 1.6;
  MeshParameters mesh_parameters{.mesh_type = "line", .lengths = {length}, .num_elements = {20}, .periodic_dims = {0}};
  mfem::Mesh mesh = buildMesh(mesh_parameters);

  constexpr int dg_order = 0;
  const int num_equations = euler::ConservativeVariables::NUM_VARS;
  Discretization dg_discretization(&mesh, dg_order, FETypes::DG, num_equations);

  constexpr double number_density = 2.1e20;
  constexpr double temperature = 281.;
  const mfem::Vector bulk_velocity({2.3, 0., 0.});
  SourceStateParameters state{.number_density = number_density, .bulk_velocity = bulk_velocity, .temperature = temperature};
  std::vector<std::unique_ptr<SourceParameters>> list_of_ic_parameters;
  list_of_ic_parameters.push_back(std::make_unique<ConstantSourceParameters>(electron_species, state));

  const LowFidelityState initial_state = buildEulerState(dg_discretization, list_of_ic_parameters);

  constexpr int es_order = 1;
  Discretization es_discretization(&mesh, es_order, FETypes::HGRAD);
  ElectrostaticFieldState es_field_state_zero(es_discretization); 

  const std::vector<Species> species_list = initial_state.getSpeciesList();
  std::vector<std::vector<std::unique_ptr<DGBC>>> empty_bcs(species_list.size());
  std::unique_ptr<LowFidelityOperations> dg_euler_operations = buildDGEulerOperations(
    dg_discretization, es_discretization, species_list, empty_bcs, empty_sources);

  const LowFidelityState rhs = dg_euler_operations->computeRHS(initial_state, es_field_state_zero);
  const LowFidelitySpeciesState& rhs_species_state = rhs.getSpeciesState(0);
  const mfem::GridFunction& rhs_grid_function = rhs_species_state.getGridFunction();
  EXPECT_EQ(0., rhs_grid_function.Norml2());
}

TEST(DGEulerOperations, computeRHS_ZeroFieldsGivesMoveIncrement) {
  constexpr double length = 0.9;
  MeshParameters mesh_parameters{.mesh_type = "line", .lengths = {length}, .num_elements = {20}, .periodic_dims = {0}};
  mfem::Mesh mesh = buildMesh(mesh_parameters);

  constexpr int dg_order = 0;
  const int num_equations = euler::ConservativeVariables::NUM_VARS;
  Discretization dg_discretization(&mesh, dg_order, FETypes::DG, num_equations);

  const mfem::Vector center{0.5 * length};
  constexpr double standard_deviation = 0.1 * length;
  const SourceStateParameters offsets{.number_density = 1.8e21, .bulk_velocity = mfem::Vector({10.2, 0., 0.})};
  const SourceStateParameters heights{.number_density = 9.2e19, .bulk_velocity = mfem::Vector({4.8, 0., 0.})};
  constexpr double pressure_offset = 29713.;
  constexpr double pressure_height = 3479.;
  std::vector<std::unique_ptr<SourceParameters>> list_of_ic_parameters;
  list_of_ic_parameters.push_back(
    std::make_unique<GaussianSourceParameters>(
      electron_species, center, standard_deviation, offsets, heights, pressure_offset, pressure_height));

  const LowFidelityState initial_state = buildEulerState(dg_discretization, list_of_ic_parameters);

  constexpr int es_order = 1;
  Discretization es_discretization(&mesh, es_order, FETypes::HGRAD);
  ElectrostaticFieldState es_field_state_zero(es_discretization); 

  const std::vector<Species> species_list = initial_state.getSpeciesList();
  std::vector<std::vector<std::unique_ptr<DGBC>>> empty_bcs(species_list.size());
  std::unique_ptr<LowFidelityOperations> dg_euler_operations = buildDGEulerOperations(
    dg_discretization, es_discretization, species_list, empty_bcs, empty_sources);

  const LowFidelityState rhs = dg_euler_operations->computeRHS(initial_state, es_field_state_zero);

  constexpr double dt = 1.5e-8;
  const LowFidelityState moved_state = dg_euler_operations->move(dt, initial_state);

  LowFidelityState updated_state(initial_state);
  updated_state.addScaledState(dt, rhs);

  for (const auto& [moved_species_state, updated_species_state] : std::views::zip(moved_state, updated_state)) {
    const mfem::GridFunction& moved_grid_function = moved_species_state.getGridFunction();
    const mfem::GridFunction& updated_grid_function = updated_species_state.getGridFunction();

    for (int i = 0; i < moved_grid_function.Size(); ++i) {
      EXPECT_DOUBLE_EQ(moved_grid_function[i], updated_grid_function[i]);
    }
  }
}

TEST(DGEulerOperations, evaluateParticleDistributionFunctionCorrectIn1D) {
  constexpr double number_density = 1e22;
  constexpr double temperature = 300;

  mfem::Vector bulk_velocity({1.0,0.0,0.0});
  const int num_elems = 1;
  constexpr int dg_order = 0;
  constexpr int num_equations = 5;
  std::shared_ptr<mfem::Mesh> mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elems));
  Discretization dg_discretization(mesh.get(), dg_order, FETypes::DG, num_equations);
  Discretization charge_discretization(mesh.get(), dg_order, FETypes::DG, num_equations);

  mfem::FiniteElementSpace finite_element_space = dg_discretization.getFeSpace();
  std::shared_ptr<DGEulerAssembly> operator_ptr = std::make_shared<DGEulerAssembly>(finite_element_space, electron_species);
  std::vector<std::shared_ptr<DGEulerAssembly>> dg_operators({operator_ptr});
  DGEulerOperations dg_euler_operations(charge_discretization, dg_operators);

  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters;
  list_of_parameters.push_back(std::make_unique<ConstantSourceParameters>(electron_species, number_density, temperature, bulk_velocity));
  LowFidelityState low_fidelity_state = buildEulerState(dg_discretization, list_of_parameters);

  const mfem::Vector position({0.5,0.0,0.0});
  const mfem::Vector particle_position(position.GetData(), 1); 
  const mfem::Vector velocity({bulk_velocity(0),0.0,0.0});
  const int low_fidelity_species_index = low_fidelity_state.getSpeciesIndex(electron_species);
  double particle_distribution_value = dg_euler_operations.evaluateParticleDistributionFunction(low_fidelity_state,particle_position,velocity,0,low_fidelity_species_index);

  mfem::Vector prim = euler::constructPrimitiveState(number_density, bulk_velocity, temperature);
  const double sigma = std::sqrt(constants::boltzmann_constant * temperature / electron_species.mass);

  const double expected_at_mean =
      number_density / std::pow(std::sqrt(2.0 * M_PI) * sigma,3);
  EXPECT_NEAR(particle_distribution_value, expected_at_mean, expected_at_mean * 1e-12);

}

TEST(DGEulerOperations, computeTotalEnergy_dg_order_0) {
  constexpr int num_elems = 10;
  constexpr double length = 2.0;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems, length);

  constexpr int dg_order = 0;
  constexpr int num_equations = euler::ConservativeVariables::NUM_VARS;
  Discretization vector_dg_discretization(&mesh, dg_order, FETypes::DG, num_equations);

  constexpr double number_density_e = 2.e17;
  const mfem::Vector bulk_velocity_e{1.2, 3.4, 5.6};
  constexpr double temperature_e = 287;

  constexpr double number_density_p = 3.e18;
  const mfem::Vector bulk_velocity_p{7.8, 9.1, 1.9};
  constexpr double temperature_p = 320;

  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters;
  list_of_parameters.push_back(
    std::make_unique<ConstantSourceParameters>(electron_species, number_density_e, temperature_e, bulk_velocity_e));
  list_of_parameters.push_back(
    std::make_unique<ConstantSourceParameters>(proton_species, number_density_p, temperature_p, bulk_velocity_p));

  LowFidelityState low_fidelity_state = buildEulerState(vector_dg_discretization, list_of_parameters);

  constexpr int es_order = 1;
  Discretization charge_discretization(&mesh, es_order, FETypes::HGRAD);
  const std::vector<Species> species_list = low_fidelity_state.getSpeciesList();
  std::vector<std::vector<std::unique_ptr<DGBC>>> empty_bcs(species_list.size());
  std::unique_ptr<LowFidelityOperations> dg_euler_operations = buildDGEulerOperations(
    vector_dg_discretization, charge_discretization, species_list, empty_bcs, empty_sources);

  const mfem::Vector primitive_state_e = euler::constructPrimitiveState(number_density_e, bulk_velocity_e, temperature_e);
  const mfem::Vector conservative_state_e = euler::convertFromPrimitiveToConservative(primitive_state_e, electron_species);
  const double total_energy_density_e = conservative_state_e[euler::ConservativeVariables::TOTAL_ENERGY_DENSITY];

  const mfem::Vector primitive_state_p = euler::constructPrimitiveState(number_density_p, bulk_velocity_p, temperature_p);
  const mfem::Vector conservative_state_p = euler::convertFromPrimitiveToConservative(primitive_state_p, proton_species);
  const double total_energy_density_p = conservative_state_p[euler::ConservativeVariables::TOTAL_ENERGY_DENSITY];

  const double expected_total_energy = (total_energy_density_e + total_energy_density_p) * length;

  const double total_energy = dg_euler_operations->computeTotalEnergy(low_fidelity_state);
  EXPECT_DOUBLE_EQ(expected_total_energy, total_energy);
}

TEST(DGEulerOperations, computeTotalKineticEnergy_ZeroVelocityGivesZeroEnergy) {
  constexpr int num_elems = 10;
  constexpr double length = 2.0;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems, length);

  constexpr int dg_order = 0;
  constexpr int num_equations = euler::ConservativeVariables::NUM_VARS;
  Discretization vector_dg_discretization(&mesh, dg_order, FETypes::DG, num_equations);

  constexpr double number_density_e = 3.e18;
  const mfem::Vector bulk_velocity_e{0., 0., 0.};
  constexpr double temperature_e = 301;

  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters;
  list_of_parameters.push_back(
    std::make_unique<ConstantSourceParameters>(electron_species, number_density_e, temperature_e, bulk_velocity_e));

  LowFidelityState low_fidelity_state = buildEulerState(vector_dg_discretization, list_of_parameters);

  constexpr int es_order = 1;
  Discretization charge_discretization(&mesh, es_order, FETypes::HGRAD);
  const std::vector<Species> species_list = low_fidelity_state.getSpeciesList();
  std::vector<std::vector<std::unique_ptr<DGBC>>> empty_bcs(species_list.size());
  std::unique_ptr<LowFidelityOperations> dg_euler_operations = buildDGEulerOperations(
    vector_dg_discretization, charge_discretization, species_list, empty_bcs, empty_sources);

  const double total_kinetic_energy = dg_euler_operations->computeTotalKineticEnergy(low_fidelity_state);
  constexpr double expected_total_kinetic_energy = 0.;
  EXPECT_DOUBLE_EQ(expected_total_kinetic_energy, total_kinetic_energy);
}

TEST(DGEulerOperations, computeTotalKineticEnergy_CorrectNonZeroKineticEnergyComputed) {
  constexpr int num_elems = 10;
  constexpr double length = 2.0;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems, length);

  constexpr int dg_order = 0;
  constexpr int num_equations = euler::ConservativeVariables::NUM_VARS;
  Discretization vector_dg_discretization(&mesh, dg_order, FETypes::DG, num_equations);

  constexpr double number_density_e = 3.e23;
  const mfem::Vector bulk_velocity_e{90.2, 170.1, 315.9};
  constexpr double temperature_e = 301;

  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters;
  list_of_parameters.push_back(
    std::make_unique<ConstantSourceParameters>(electron_species, number_density_e, temperature_e, bulk_velocity_e));

  LowFidelityState low_fidelity_state = buildEulerState(vector_dg_discretization, list_of_parameters);

  constexpr int es_order = 1;
  Discretization charge_discretization(&mesh, es_order, FETypes::HGRAD);
  const std::vector<Species> species_list = low_fidelity_state.getSpeciesList();
  std::vector<std::vector<std::unique_ptr<DGBC>>> empty_bcs(species_list.size());
  std::unique_ptr<LowFidelityOperations> dg_euler_operations = buildDGEulerOperations(
    vector_dg_discretization, charge_discretization, species_list, empty_bcs, empty_sources);

  const double total_kinetic_energy = dg_euler_operations->computeTotalKineticEnergy(low_fidelity_state);

  const mfem::Vector primitive_state_e = euler::constructPrimitiveState(number_density_e, bulk_velocity_e, temperature_e);
  const mfem::Vector conservative_state_e = euler::convertFromPrimitiveToConservative(primitive_state_e, electron_species);
  const double kinetic_energy_density = euler::getKineticEnergyDensityFromConservativeState(conservative_state_e);
  const double expected_total_kinetic_energy = kinetic_energy_density * length;
  EXPECT_DOUBLE_EQ(expected_total_kinetic_energy, total_kinetic_energy);
}

TEST(DGEulerOperations, LowFidelityMomentsAreCorrectForMaxwellianIn3D) {
  Species species{.charge = -constants::elementary_charge, .mass = constants::electron_mass};
  constexpr double number_density = 1e22;
  constexpr double temperature = 300;
  mfem::Vector bulk_velocity({293.0,581.0,902.0});

  const int num_elems = 5;
  constexpr int dg_order = 0;
  constexpr int num_equations = 5;
  constexpr mfem::Element::Type element_type = mfem::Element::HEXAHEDRON;
  std::shared_ptr<mfem::Mesh> mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(
    num_elems,
    num_elems,
    num_elems,
    element_type
  ));

  Discretization dg_discretization(mesh.get(), dg_order, FETypes::DG, num_equations);

  constexpr int charge_order = 1;
  Discretization charge_discretization(mesh.get(), charge_order, FETypes::HGRAD);

  mfem::FiniteElementSpace finite_element_space = dg_discretization.getFeSpace();
  std::shared_ptr<DGEulerAssembly> operator_ptr = std::make_shared<DGEulerAssembly>(finite_element_space, species);
  std::vector<std::shared_ptr<DGEulerAssembly>> dg_operators({operator_ptr});
  DGEulerOperations dg_euler_operations(charge_discretization, dg_operators);

  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters;
  list_of_parameters.push_back(std::make_unique<ConstantSourceParameters>(species, number_density, temperature,bulk_velocity));
  LowFidelityState low_fidelity_state = buildEulerState(dg_discretization, list_of_parameters);

  std::unordered_map<Species, mfem::Vector> computed_number_density = dg_euler_operations.getCellAveragedNumberDensity(low_fidelity_state);
  std::unordered_map<Species, mfem::DenseMatrix> computed_bulk_velocity = dg_euler_operations.getCellAveragedBulkVelocity(low_fidelity_state);
  std::unordered_map<Species, mfem::Vector> computed_temperature = dg_euler_operations.getCellAveragedTemperature(low_fidelity_state);

  for (int elem_id = 0; elem_id < num_elems; ++elem_id)
  {
    EXPECT_DOUBLE_EQ(computed_number_density.at(species)(elem_id),number_density);
    EXPECT_DOUBLE_EQ(computed_bulk_velocity.at(species)(0,elem_id),bulk_velocity(0));
    EXPECT_DOUBLE_EQ(computed_bulk_velocity.at(species)(1,elem_id),bulk_velocity(1));
    EXPECT_DOUBLE_EQ(computed_bulk_velocity.at(species)(2,elem_id),bulk_velocity(2));
    EXPECT_DOUBLE_EQ(computed_temperature.at(species)(elem_id),temperature);
  }
}

} // namespace
