#include <libmfpic/Constants.hpp>
#include <libmfpic/DGGhostBoundaryIntegrator.hpp>
#include <libmfpic/DGEulerBoundaryConditions.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/KineticFluxBC.hpp>

#include <gtest/gtest.h>
#include <libmfpic/Euler.hpp>
#include <libmfpic/KineticFluxBC.hpp>
#include <mfem/fem/eltrans.hpp>
#include <mfem/fem/hyperbolic.hpp>
#include <mfem/fem/lininteg.hpp>
#include <mfem/linalg/densemat.hpp>
#include <random>

namespace {

using namespace mfpic;

TEST(DGEulerBoundaryConditions, DGEulerReflectingBCSetsGhostCorrectly) {
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(10);
  DGEulerReflectingBC bc;

  constexpr int num_eqns = 5;
  constexpr int num_dof  = 4;

  mfem::Vector normal {.65, -.2, .11};
  normal /= normal.Norml2();

  mfem::DenseMatrix in(num_dof, num_eqns), out(num_dof,num_eqns);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> pos_dis(.1, 2.0);
  std::uniform_real_distribution<> dis(-2.0, 2.0);

  for (int idof = 0; idof < num_dof; ++idof) {
    in(idof, 0) = pos_dis(gen);
    in(idof, 1) = dis(gen);
    in(idof, 2) = dis(gen);
    in(idof, 3) = dis(gen);
    in(idof, 4) = pos_dis(gen);
  }

  bc.setDOFsInGhost(in, normal, out);

  for (int idof = 0; idof < num_dof; ++idof) {
    const mfem::Vector momentum {in(idof,1), in(idof,2), in(idof,3)};
    mfem::Vector normal_momentum = normal;
    normal_momentum *= (momentum * normal);
    EXPECT_DOUBLE_EQ(out(idof,0), in(idof,0));
    EXPECT_DOUBLE_EQ(out(idof,4), in(idof,4));
    EXPECT_DOUBLE_EQ(out(idof,1), in(idof,1) - 2. * normal_momentum(0));
    EXPECT_DOUBLE_EQ(out(idof,2), in(idof,2) - 2. * normal_momentum(1));
    EXPECT_DOUBLE_EQ(out(idof,3), in(idof,3) - 2. * normal_momentum(2));
  }
}

/// numerical flux that returns F(U) \cdot n for 1 or 2
class PickOneFlux : public mfem::NumericalFlux {
  public:
    PickOneFlux(const mfem::FluxFunction &flux_function, const bool pick_ghost) :
      mfem::NumericalFlux(flux_function),
      pick_ghost_(pick_ghost) {};

    mfem::real_t Eval(const mfem::Vector &state1, const mfem::Vector &state2,
                      const mfem::Vector &nor, mfem::FaceElementTransformations &transformations,
                      mfem::Vector &flux) const override 
    {
      const mfem::real_t speed = pick_ghost_ ? 
        fluxFunction.ComputeFluxDotN(state2, nor, transformations, flux) :
        fluxFunction.ComputeFluxDotN(state1, nor, transformations, flux);
      return speed;
    };
  private:
    const bool pick_ghost_;
};
/// linear flux for testing, only has a boundary term 
class LinearFlux : public mfem::FluxFunction {
  public:
    LinearFlux(int dim) :
      mfem::FluxFunction(5, dim) {}
    mfem::real_t ComputeFluxDotN(const mfem::Vector & state, const mfem::Vector & normal,
                                 mfem::FaceElementTransformations &, 
                                 mfem::Vector & flux_dot_n) const override 
    {
      const double darea = normal.Norml2(); // mfem uses the normal to get the weight correct 
      mfem::Vector unit_normal = normal;
      unit_normal /= normal.Norml2();
      const mfem::Vector momentum(state.GetData() + 1, dim);
      const mfem::real_t total_energy = state(1 + dim);
      mfem::Vector normal_momentum = unit_normal;
      normal_momentum *= (momentum * unit_normal);
      flux_dot_n(0) = momentum * unit_normal;
      for (int d = 0; d < dim; ++d)
        flux_dot_n(1 + d) = normal_momentum(d);
      flux_dot_n(1 + dim) = total_energy;
      flux_dot_n *= darea;
      return momentum * unit_normal * darea;
    }
    mfem::real_t ComputeFlux(const mfem::Vector & , mfem::ElementTransformation & ,
                             mfem::DenseMatrix & ) const override {return 0;}
};

TEST(DGEulerBoundaryConditions, DGEulerReflectingBCCheckGhostBoundaryIntegrator) {
  constexpr double tolerance = 1e-12;
  int order = 1;
  int dim = 3; 
  int num_equations = 5;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(3, 3, 3, mfem::Element::HEXAHEDRON, 1, 1, 1);
  Discretization discretization(&mesh, order, FETypes::DG, num_equations);
  mfem::FiniteElementSpace finite_element_space = discretization.getFeSpace();

  constexpr mfem::real_t c0(12.7), c1(-9.4), c2(2.2), c3(9.1);
  auto solution_vec = [&](const mfem::Vector &x, mfem::Vector &y) { 
    mfem::real_t base_val = c0 + c1 * x[0] + c2 * x[1] + c3 * x[2];
    for (int i = 0; i < 5; ++i)
      y[i] = base_val * (i + 1);
  };
  mfem::VectorFunctionCoefficient fluid_coeff(num_equations,solution_vec);
  mfem::GridFunction fluid_dofs(&discretization.getFeSpace());
  fluid_dofs.ProjectCoefficient(fluid_coeff);

  LinearFlux linear_flux(dim);
  constexpr bool pick_ghost = true;
  PickOneFlux numerical_flux(linear_flux, pick_ghost);
  mfem::NonlinearForm form(&finite_element_space);
  constexpr int boundary_attribute = 3; // this should be the right face (x = 1)

<<<<<<< HEAD
DGGhostBC bc(boundary_attribute, mesh, Species(), std::make_unique<DGEulerReflectingBC>());
=======
  DGGhostBC bc(boundary_attribute, mesh, Species(), std::make_unique<DGEulerReflectingBC>());
>>>>>>> origin/main
  auto integrator = bc.makeIntegrator(numerical_flux);
  form.AddBdrFaceIntegrator(integrator.release(), bc.boundary_attribute_has_boundary_condition);

  mfem::Vector rhs(fluid_dofs.Size());
  rhs = 0.;
  form.Mult(fluid_dofs, rhs);

  auto f_dot_n_expected = [&](const mfem::Vector &x, mfem::Vector &y) { 
    mfem::Vector u(5);
    solution_vec(x,u);
    mfem::Vector normal{1.,0.,0.};

    // first reverse the momentum
    {
      const double p_dot_n = normal[0] * u[1] + normal[1] * u[2] + normal[2] * u[3];

      mfem::Vector normal_momentum = normal;
      normal_momentum *= p_dot_n;
    
      u[1] -= 2 * normal_momentum[0];
      u[2] -= 2 * normal_momentum[1];
      u[3] -= 2 * normal_momentum[2];
    }
    y = 0.;
    // now compute - F(U) dot n (negative due to contribution to left cell)
    {
      const double p_dot_n = normal[0] * u[1] + normal[1] * u[2] + normal[2] * u[3];

      mfem::Vector normal_momentum = normal;
      normal_momentum *= p_dot_n;
      y[0] = -p_dot_n;
      y[1] = -normal_momentum[0];
      y[2] = -normal_momentum[1];
      y[3] = -normal_momentum[2];
      y[4] = -u[4];
    }
  };

  mfem::VectorFunctionCoefficient exact_coeff(num_equations,f_dot_n_expected);
  mfem::LinearForm exact_form(&discretization.getFeSpace());
  exact_form.AddBdrFaceIntegrator(
    new mfem::VectorBoundaryLFIntegrator(exact_coeff), bc.boundary_attribute_has_boundary_condition);
  exact_form.Assemble();

  mfem::Vector diff = exact_form;
  auto exact_l2 = diff.Norml2();
  diff.Add(-1.0, rhs);
  auto diff_l2 = diff.Norml2();

  EXPECT_NEAR(diff_l2/exact_l2, 0., tolerance);

}

constexpr Species electron_species{.charge = -constants::elementary_charge, .mass = constants::electron_mass, .specific_heat_ratio = 5./3.};

TEST(DGEulerBoundaryConditions, DummyFluxFunctionHasNoComputeFlux) {
  auto dummy_flux = DummyFluxFunction(5, 3);
  auto dummy_transform = mfem::IsoparametricTransformation();
  auto dummy_matrix = mfem::DenseMatrix(); 
  EXPECT_DEATH(dummy_flux.ComputeFlux(mfem::Vector(), dummy_transform, dummy_matrix),
               "ComputeFlux cannot be called!");
}

mfem::Vector expectedFDotN(const mfem::Vector conservative_state, const mfem::Vector unit_normal, const Species species)
{
  // get orthonormal basis
  mfem::Vector e(3);
  if (std::abs(unit_normal(0)) < .9) {
    e = mfem::Vector({1.0,0.0,0.0});
  } else {
    e = mfem::Vector({0.0,1.0,0.0});
  }
  mfem::Vector tangent_1(3), tangent_2(3);
  unit_normal.cross3D(e, tangent_1);
  tangent_1 /= tangent_1.Norml2();
  unit_normal.cross3D(tangent_1,tangent_2);
  tangent_2 /= tangent_2.Norml2();

  mfem::Vector f_dot_n(5);
  const double temperature = euler::getTemperatureFromConservativeState(conservative_state, species);
  const double number_density = euler::getNumberDensityFromConservativeState(conservative_state, species);

  const double u_m = std::sqrt(2. * constants::boltzmann_constant * temperature / species.mass);
  mfem::Vector reduced_velocity = euler::getBulkVelocityFromConservativeState(conservative_state);
  reduced_velocity /= u_m;

  const double reduced_dot_normal = reduced_velocity * unit_normal;
  const double reduced_dot_t1     = reduced_velocity * tangent_1;
  const double reduced_dot_t2     = reduced_velocity * tangent_2;

  // see sec 1.3.6 in Boyd and Schwartzentruber
  const double density_flux = 1. / 4. * conservative_state[euler::ConservativeVariables::MASS_DENSITY] * 
    std::sqrt(8. * constants::boltzmann_constant * temperature  / (M_PI * species.mass)) * (
      std::exp(-reduced_dot_normal*reduced_dot_normal) + 
      std::sqrt(M_PI) * reduced_dot_normal * (std::erf(reduced_dot_normal) + 1.)
    );
  const double normal_momentum_flux = conservative_state[euler::ConservativeVariables::MASS_DENSITY] * 
    constants::boltzmann_constant * temperature / species.mass * (
      reduced_dot_normal / std::sqrt(M_PI) * std::exp(-reduced_dot_normal*reduced_dot_normal) +
      (.5 + reduced_dot_normal * reduced_dot_normal) * (std::erf(reduced_dot_normal) + 1.) 
    );
  const double t1_momentum_flux = density_flux * reduced_dot_t1;
  const double t2_momentum_flux = density_flux * reduced_dot_t2;
  const double energy_flux = 1. / 4. * number_density *  constants::boltzmann_constant * temperature * 
  std::sqrt(8. * constants::boltzmann_constant * temperature  / (M_PI * species.mass)) * (
    (reduced_velocity * reduced_velocity + 2.) * std::exp(-reduced_dot_normal*reduced_dot_normal) +
    (reduced_velocity * reduced_velocity + 5./2.) * std::sqrt(M_PI) * reduced_dot_normal * (std::erf(reduced_dot_normal) + 1.)
  );

  f_dot_n(0) = density_flux;
  for (int d = 0; d < 3; ++d)
    f_dot_n(1 + d) = normal_momentum_flux * unit_normal(d) + t1_momentum_flux * tangent_1(d) + t2_momentum_flux * tangent_2(d);
  f_dot_n(4) = energy_flux;

  return f_dot_n;
}

TEST(DGEulerBoundaryConditions, KineticFluxFluxFunctionComputesCorrectFluxIn1D) {

  constexpr double tol = 1e-12;
  constexpr int spatial_dim = 1;
  auto dummy_flux = DummyFluxFunction(5, spatial_dim);
  KineticFluxNumericalFlux flux(dummy_flux, electron_species);
  auto dummy_transform = mfem::FaceElementTransformations();

  const mfem::Vector normal {2.567, 0., 0.}; // mfem's normals are not necessarily unit vectors
  mfem::Vector unit_normal = normal;
  unit_normal /= normal.Norml2();

  constexpr double temperature = 293.;
  const mfem::Vector bulk_velocity {12.345, 0., 0.};
  constexpr double number_density = 1.567e16;

  mfem::Vector primitive_state = euler::constructPrimitiveState(number_density, bulk_velocity, temperature);
  mfem::Vector conservative_state = euler::convertFromPrimitiveToConservative(primitive_state, electron_species);
  const mfem::Vector f_dot_n_expected = expectedFDotN(conservative_state, unit_normal, electron_species);
  mfem::Vector f_dot_n(5);

  flux.Eval(conservative_state, conservative_state, normal, dummy_transform, f_dot_n);

  for (int i = 0; i < 5; ++i)
    EXPECT_NEAR(f_dot_n(i), f_dot_n_expected(i) * normal.Norml2(), tol);

}

TEST(DGEulerBoundaryConditions, KineticFluxFluxFunctionComputesCorrectFluxIn2D) {

  constexpr double tol = 1e-7;
  constexpr int spatial_dim = 2;
  auto dummy_flux = DummyFluxFunction(5, spatial_dim);
  KineticFluxNumericalFlux flux(dummy_flux, electron_species);
  auto dummy_transform = mfem::FaceElementTransformations();

  const mfem::Vector normal {2.567, -1.234, 0.}; // mfem's normals are not necessarily unit vectors
  mfem::Vector unit_normal = normal;
  unit_normal /= normal.Norml2();

  constexpr double temperature = 293.;
  const mfem::Vector bulk_velocity {12.345, -8.924, 0.};
  constexpr double number_density = 1.567e16;

  mfem::Vector primitive_state = euler::constructPrimitiveState(number_density, bulk_velocity, temperature);
  mfem::Vector conservative_state = euler::convertFromPrimitiveToConservative(primitive_state, electron_species);
  const mfem::Vector f_dot_n_expected = expectedFDotN(conservative_state, unit_normal, electron_species);
  mfem::Vector f_dot_n(5);

  flux.Eval(conservative_state, conservative_state, normal, dummy_transform, f_dot_n);

  for (int i = 0; i < 5; ++i)
    EXPECT_NEAR(f_dot_n(i), f_dot_n_expected(i) * normal.Norml2(), tol);

}

TEST(DGEulerBoundaryConditions, KineticFluxFluxFunctionComputesCorrectFluxIn3D) {

  constexpr double tol = 1e-7;
  constexpr int spatial_dim = 3;
  auto dummy_flux = DummyFluxFunction(5, spatial_dim);
  KineticFluxNumericalFlux flux(dummy_flux, electron_species);
  auto dummy_transform = mfem::FaceElementTransformations();

  const mfem::Vector normal {2.567, -1.234, 3.571}; // mfem's normals are not necessarily unit vectors
  mfem::Vector unit_normal = normal;
  unit_normal /= normal.Norml2();

  constexpr double temperature = 293.;
  const mfem::Vector bulk_velocity {12.345, -8.924, -17.234};
  constexpr double number_density = 1.567e16;

  mfem::Vector primitive_state = euler::constructPrimitiveState(number_density, bulk_velocity, temperature);
  mfem::Vector conservative_state = euler::convertFromPrimitiveToConservative(primitive_state, electron_species);
  const mfem::Vector f_dot_n_expected = expectedFDotN(conservative_state, unit_normal, electron_species);
  mfem::Vector f_dot_n(5);

  flux.Eval(conservative_state, conservative_state, normal, dummy_transform, f_dot_n);

  for (int i = 0; i < 5; ++i)
    EXPECT_NEAR(f_dot_n(i), f_dot_n_expected(i) * normal.Norml2(), tol) << i;

}

TEST(DGEulerBoundaryConditions, KineticFluxFluxFunctionCorrectHighMachBehavior) {

  constexpr double tol = 1e-12;
  constexpr int spatial_dim = 1;
  auto dummy_flux = DummyFluxFunction(5, spatial_dim);
  KineticFluxNumericalFlux flux(dummy_flux, electron_species);  
  auto dummy_transform = mfem::FaceElementTransformations();

  const mfem::Vector normal {2.567, 0., 0.}; // mfem's normals are not necessarily unit vectors
  mfem::Vector unit_normal = normal;
  unit_normal /= normal.Norml2();

  constexpr double mach = 10.;

  const mfem::Vector bulk_velocity {10000.0, 0., 0.};
  const double temperature = bulk_velocity(0) * bulk_velocity(0) * electron_species.mass / (mach * mach * electron_species.specific_heat_ratio * constants::boltzmann_constant);
  constexpr double number_density = 1.567e16;

  mfem::Vector primitive_state = euler::constructPrimitiveState(number_density, bulk_velocity, temperature);
  mfem::Vector conservative_state = euler::convertFromPrimitiveToConservative(primitive_state, electron_species);
  const double pressure = euler::getPressureFromConservativeState(conservative_state, electron_species);
  mfem::Vector f_dot_n_expected(5);
  // should just get right-going flux
  f_dot_n_expected(0) = conservative_state[euler::ConservativeVariables::X_MOMENTUM_DENSITY];
  f_dot_n_expected(1) = f_dot_n_expected(0) * primitive_state[euler::PrimitiveVariables::X_BULK_VELOCITY] + pressure;
  f_dot_n_expected(2) = 0.;
  f_dot_n_expected(3) = 0.;
  f_dot_n_expected(4) = primitive_state[euler::PrimitiveVariables::X_BULK_VELOCITY] * (pressure + conservative_state[euler::ConservativeVariables::TOTAL_ENERGY_DENSITY]);
  mfem::Vector f_dot_n(5);

  flux.Eval(conservative_state, conservative_state, normal, dummy_transform, f_dot_n);

  for (int i = 0; i < 5; ++i)
    EXPECT_NEAR(f_dot_n(i), f_dot_n_expected(i) * normal.Norml2(), tol);
}

TEST(DGEulerBoundaryConditions, KineticFluxFluxFunctionCorrectLeftGoingExtreme) {

  constexpr double tol = 1e-12;
  constexpr int spatial_dim = 1;
  auto dummy_flux = DummyFluxFunction(5, spatial_dim);
  KineticFluxNumericalFlux flux(dummy_flux, electron_species);
  auto dummy_transform = mfem::FaceElementTransformations();

  const mfem::Vector normal {2.567, 0., 0.}; // mfem's normals are not necessarily unit vectors
  mfem::Vector unit_normal = normal;
  unit_normal /= normal.Norml2();

  constexpr double mach = 10.;

  const mfem::Vector bulk_velocity {-10000.0, 0., 0.};
  const double temperature = bulk_velocity(0) * bulk_velocity(0) * electron_species.mass / (mach * mach * electron_species.specific_heat_ratio * constants::boltzmann_constant);
  constexpr double number_density = 1.567e16;

  mfem::Vector primitive_state = euler::constructPrimitiveState(number_density, bulk_velocity, temperature);
  mfem::Vector conservative_state = euler::convertFromPrimitiveToConservative(primitive_state, electron_species);
  mfem::Vector f_dot_n_expected(5);
  // should get no flux
  f_dot_n_expected = 0.;
  mfem::Vector f_dot_n(5);

  flux.Eval(conservative_state, conservative_state, normal, dummy_transform, f_dot_n);

  for (int i = 0; i < 5; ++i)
    EXPECT_NEAR(f_dot_n(i), f_dot_n_expected(i) * normal.Norml2(), tol);
}

TEST(DGEulerBoundaryConditions, KineticFluxFluxFunctionCorrectNoBulk) {

  constexpr double tol = 1e-12;
  constexpr int spatial_dim = 1;
  auto dummy_flux = DummyFluxFunction(5, spatial_dim);
  KineticFluxNumericalFlux flux(dummy_flux, electron_species);
  auto dummy_transform = mfem::FaceElementTransformations();

  const mfem::Vector normal {2.567, 0., 0.}; // mfem's normals are not necessarily unit vectors
  mfem::Vector unit_normal = normal;
  unit_normal /= normal.Norml2();

  constexpr double temperature = 293.;
  const mfem::Vector bulk_velocity {0., 0., 0.};
  constexpr double number_density = 1.567e16;

  const double thermal_speed = std::sqrt(2. * constants::boltzmann_constant * temperature / (M_PI * electron_species.mass));

  mfem::Vector primitive_state = euler::constructPrimitiveState(number_density, bulk_velocity, temperature);
  mfem::Vector conservative_state = euler::convertFromPrimitiveToConservative(primitive_state, electron_species);
  const double pressure = euler::getPressureFromConservativeState(conservative_state, electron_species);
  mfem::Vector f_dot_n_expected(5);
  // should just get thermal fluctuations
  // see sec 1.3.6 in Boyd and Schwartzentruber
  f_dot_n_expected(0) = conservative_state[euler::ConservativeVariables::MASS_DENSITY] * thermal_speed / 2.;
  f_dot_n_expected(1) = pressure / 2.;
  f_dot_n_expected(2) = 0.;
  f_dot_n_expected(3) = 0.;
  f_dot_n_expected(4) = thermal_speed * pressure; 
  mfem::Vector f_dot_n(5);

  flux.Eval(conservative_state, conservative_state, normal, dummy_transform, f_dot_n);

  for (int i = 0; i < 5; ++i)
    EXPECT_NEAR(f_dot_n(i), f_dot_n_expected(i) * normal.Norml2(), tol) << i;
}

} // namespace
