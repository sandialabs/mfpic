#include <libmfpic/ConjugateGradientLinearSolver.hpp>
#include <libmfpic/Constants.hpp>
#include <libmfpic/DirichletBoundaryConditions.hpp>
#include <libmfpic/DirichletBoundaryConditionsFactory.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/ElectrostaticFieldState.hpp>
#include <libmfpic/LoadUniformMaxwellianParticles.hpp>
#include <libmfpic/ParticleOperations.hpp>
#include <libmfpic/ReflectingParticleBoundary.hpp>

#include <mfem/mfem.hpp>

#include <gtest/gtest.h>

#include <random>

namespace {

using namespace mfpic;

TEST(ConjugateGradientLinearSolver, SolvingIdentityGivesRHS) {
  constexpr int size = 10;
  mfem::Vector identity_diagonal(size);
  identity_diagonal = 1.;
  mfem::SparseMatrix identity_matrix(identity_diagonal);

  ConjugateGradientLinearSolver cg_linear_solver;

  mfem::Vector rhs_vector(size);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(1.0, 2.0);
  for (int i = 0; i < size; ++i) {
    rhs_vector[i] = dis(gen);
  }

  mfem::Vector sol_vector(size);
  cg_linear_solver.solve(identity_matrix, sol_vector, rhs_vector);

  for (int i = 0; i < size; ++i) {
    EXPECT_NEAR(sol_vector[i], rhs_vector[i], 1e-12);
  }
}

double computeResidualNorm(const mfem::SparseMatrix& matrix, const mfem::Vector& sol, const mfem::Vector& rhs) {
  mfem::Vector residual(sol.Size());
  matrix.Mult(sol, residual);
  residual -= rhs;
  const double residual_norm = residual.Norml2();
  return residual_norm;
}

mfem::Vector generateRandomVector(const int size) {
  std::mt19937 gen;
  std::uniform_real_distribution<> dis(1.0, 10.0);

  mfem::Vector random_vector(size);
  for (int i = 0; i < size; ++i) {
    random_vector[i] = dis(gen);
  }

  return random_vector;
}

TEST(ConjugateGradientLinearSolver, InvertDiffusionMatrix) {
  constexpr int size = 10;
  constexpr int num_elems = size - 1;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);

  constexpr int hgrad_order = 1;
  mfem::H1_FECollection hgrad_finite_element_collection(hgrad_order, mesh.Dimension());
  mfem::FiniteElementSpace hgrad_finite_element_space(&mesh, &hgrad_finite_element_collection);

  mfem::BilinearForm laplacian_bilinear_form(&hgrad_finite_element_space);
  laplacian_bilinear_form.AddDomainIntegrator(new mfem::DiffusionIntegrator);
  laplacian_bilinear_form.Assemble();

  mfem::SparseMatrix laplacian_matrix;
  mfem::Array<int> dirichlet_boundary_indices{0};
  laplacian_bilinear_form.FormSystemMatrix(dirichlet_boundary_indices, laplacian_matrix);

  mfem::Vector rhs_vector = generateRandomVector(size);
  mfem::Vector sol_vector(size);
  sol_vector = 0.;

  const double initial_residual_norm = computeResidualNorm(laplacian_matrix, sol_vector, rhs_vector);

  ConjugateGradientLinearSolver cg_linear_solver;
  cg_linear_solver.solve(laplacian_matrix, sol_vector, rhs_vector);

  const double final_residual_norm = computeResidualNorm(laplacian_matrix, sol_vector, rhs_vector);
  const double relative_residual_norm = final_residual_norm / initial_residual_norm;
  EXPECT_LT(relative_residual_norm, 1e-12);
}

mfem::SparseMatrix generateRandomSPDMatrix(const int size) {
  const int row_size = size;
  mfem::SparseMatrix random_matrix(size, size, row_size);

  std::mt19937 gen;
  std::uniform_real_distribution<> dis(1.0, 10.0);
  for (int i = 0; i < size; ++i) {
    mfem::Array<int> col_indices(size);
    mfem::Vector row_entries(size);
    for (int j = 0; j < size; ++j) {
      col_indices[j] = j;
      row_entries[j] = dis(gen);
    }
    random_matrix.SetRow(i, col_indices, row_entries);
  }

  mfem::SparseMatrix random_spd_matrix(size, size, row_size);
  for (int i = 0; i < size; ++i) {
    mfem::Array<int> col_indices(size);
    mfem::Vector row_entries(size);
    for (int j = 0; j < size; ++j) {
      col_indices[j] = j;
      mfem::real_t value = 0.0;
      for (int k = 0; k < size; ++k) {
        value += random_matrix(i, k) * random_matrix(j, k);
      }
      row_entries[j] = value;
    }
    random_spd_matrix.SetRow(i, col_indices, row_entries);
  }

  return random_spd_matrix;
}

TEST(ConjugateGradientLinearSolver, InvertRandomSPDMatrix) {
  constexpr int size = 10;
  const mfem::SparseMatrix random_spd_matrix = generateRandomSPDMatrix(size);


  mfem::Vector rhs_vector = generateRandomVector(size);
  mfem::Vector sol_vector(size);
  sol_vector = 0.;

  const double initial_residual_norm = computeResidualNorm(random_spd_matrix, sol_vector, rhs_vector);

  ConjugateGradientLinearSolver cg_linear_solver;
  cg_linear_solver.solve(random_spd_matrix, sol_vector, rhs_vector);

  const double final_residual_norm = computeResidualNorm(random_spd_matrix, sol_vector, rhs_vector);
  const double relative_residual_norm = final_residual_norm / initial_residual_norm;
  EXPECT_LT(relative_residual_norm, 1e-11);
}

TEST(ConjugateGradientLinearSolver, InvertElectrostaticMatrix) {
  constexpr double number_density = 1e16;
  constexpr double qe = -constants::elementary_charge;
  constexpr double me = constants::electron_mass;
  const double plasma_frequency = sqrt((number_density * qe * qe) / (me * constants::permittivity));
  const double dt = 2. * M_PI / (120 * plasma_frequency);
  constexpr double bulk_speed = 3.2e5;
  const double dx = bulk_speed * dt;
  constexpr int num_elems = 400;
  const double length = num_elems * dx;
  auto mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elems, length));

  constexpr int electrostatic_order = 1;
  Discretization electrostatic_discretization(mesh.get(), electrostatic_order, FETypes::HGRAD);

  mfem::BilinearForm electrostatic_bilinear_form(&electrostatic_discretization.getFeSpace());
  mfem::ConstantCoefficient permittivity(constants::permittivity);
  electrostatic_bilinear_form.AddDomainIntegrator(new mfem::DiffusionIntegrator(permittivity));
  electrostatic_bilinear_form.Assemble();

  std::unique_ptr<DirichletBoundaryConditions> dirichlet_bcs = buildDirichletBoundaryConditions({}, electrostatic_discretization);

  ElectrostaticFieldState field_state(electrostatic_discretization);
  mfem::GridFunction potential = field_state.getPotential();

  auto default_particle_boundary_factory = std::make_shared<ReflectingParticleBoundaryFactory>();
  ParticleOperations particle_operations(
    electrostatic_discretization,
    {},
    default_particle_boundary_factory);

  Species electron_species{.charge = qe, .mass = me};
  Species proton_species{.charge = constants::elementary_charge, .mass = constants::proton_mass};
  constexpr int num_particles = 100000;
  constexpr double temperature = 0.;
  const mfem::Vector bulk_velocity{0, 0, 0};
  std::default_random_engine generator;
  ParticleContainer particle_container;

  particle_container.addParticles(loadUniformMaxwellianParticles(
    electron_species,
    bulk_velocity,
    temperature,
    number_density,
    num_particles,
    generator,
    mesh));
  particle_container.addParticles(loadUniformMaxwellianParticles(
    proton_species,
    bulk_velocity,
    temperature,
    number_density,
    num_particles,
    generator,
    mesh));

  IntegratedCharge integrated_charge = particle_operations.assembleCharge(particle_container);
  mfem::Vector integrated_charge_vector = integrated_charge.getIntegratedCharge();

  mfem::Array<int> dirichlet_dof_indices = dirichlet_bcs->getDirichletBoundaryDofIndices();
  mfem::Vector solution_vector;
  mfem::Vector rhs_vector;
  mfem::SparseMatrix negative_eps_laplace_matrix;
  electrostatic_bilinear_form.FormLinearSystem(
    dirichlet_dof_indices,
    potential,
    integrated_charge_vector,
    negative_eps_laplace_matrix,
    solution_vector,
    rhs_vector);

  const double initial_residual_norm = computeResidualNorm(negative_eps_laplace_matrix, solution_vector, rhs_vector);

  constexpr double relative_tolerance = 1e-8;
  constexpr double absolute_tolerance = 0;
  ConjugateGradientLinearSolver cg_linear_solver(relative_tolerance, absolute_tolerance);
  cg_linear_solver.solve(negative_eps_laplace_matrix, solution_vector, rhs_vector);

  const double final_residual_norm = computeResidualNorm(negative_eps_laplace_matrix, solution_vector, rhs_vector);
  const double norm_ratio = final_residual_norm / initial_residual_norm;
  EXPECT_LT(norm_ratio, relative_tolerance);
}

}