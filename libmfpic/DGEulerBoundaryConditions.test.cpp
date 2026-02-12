#include <libmfpic/DGGhostBoundaryIntegrator.hpp>
#include <libmfpic/DGEulerBoundaryConditions.hpp>
#include <libmfpic/Discretization.hpp>

#include <gtest/gtest.h>
#include <mfem/fem/eltrans.hpp>
#include <mfem/fem/hyperbolic.hpp>
#include <mfem/fem/lininteg.hpp>
#include <mfem/linalg/densemat.hpp>
#include <random>

namespace {

using namespace mfpic;

TEST(DGEulerBoundaryConditions, DGEulerReflectingBCSetsGhostCorrectly) {
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(10);
  DGEulerReflectingBC bc(1, mesh);

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
    EXPECT_FLOAT_EQ(out(idof,0), in(idof,0));
    EXPECT_FLOAT_EQ(out(idof,4), in(idof,4));
    EXPECT_FLOAT_EQ(out(idof,1), in(idof,1) - 2. * normal_momentum(0));
    EXPECT_FLOAT_EQ(out(idof,2), in(idof,2) - 2. * normal_momentum(1));
    EXPECT_FLOAT_EQ(out(idof,3), in(idof,3) - 2. * normal_momentum(2));
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
  DGEulerReflectingBC bc(boundary_attribute, mesh);
  form.AddBdrFaceIntegrator(new DGGhostBoundaryIntegrator(numerical_flux, bc), bc.boundary_attribute_has_boundary_condition);

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

};

} // namespace
