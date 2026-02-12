#include <libmfpic/Discretization.hpp>
#include <libmfpic/DGEulerAssembly.hpp>
#include <libmfpic/DGEulerKEIntegrator.hpp>
#include <libmfpic/DGEulerMaxwellSourceIntegrator.hpp>
#include <libmfpic/ElectrostaticFieldState.hpp>
#include <libmfpic/Species.hpp>
#include <libmfpic/LowFidelityState.hpp>
#include <gtest/gtest.h>
#include <mfem.hpp>

namespace {

using namespace mfpic;

TEST(DGEulerAssembly, EulerFluxIsCorrectForConstantStateIn3D) {
  int order = 1;
  int dim = 3;
  int num_equations = 5;
  double specific_heat_ratio = 1.4;
  const Species species{.specific_heat_ratio = specific_heat_ratio};
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(3, 3, 3, mfem::Element::HEXAHEDRON, 1, 1, 1);
  Discretization discretization(&mesh, order, FETypes::DG, num_equations);
  mfem::FiniteElementSpace finite_element_space = discretization.getFeSpace();
  DGEulerAssembly dg_euler_operator(finite_element_space, species);

  double rho = 54.0;
  double u = 328.0;
  double v = 92.0;
  double w = 31.0;
  double p = 103.0;
  mfem::Vector velocity_field{u,v,w};
  double e = p / ((species.specific_heat_ratio - 1.0)*rho) + 0.5 * (u*u + v*v + w*w);
  double enthalpy = e + p/rho;

  const mfem::Vector fluid_constant_vals{rho,rho*u,rho*v,rho*w,rho*e};
  mfem::VectorConstantCoefficient fluid_constant_coeff(fluid_constant_vals);
  mfem::GridFunction fluid_dofs(&finite_element_space);
  fluid_dofs.ProjectCoefficient(fluid_constant_coeff);

  mfem::DenseMatrix flux_in_current_element;
  mfem::DenseMatrix state_in_current_element;
  mfem::Vector state_at_current_dof(num_equations);
  mfem::Array<int> vector_dofs;
  mfem::Vector xval,tempval;
  mfem::DenseMatrix flux_for_current_equation;

  for (int element=0; element<finite_element_space.GetNE(); element++) {
    int num_dof = finite_element_space.GetFE(element)->GetDof(); // 2 DOFs per element
    finite_element_space.GetElementVDofs(element, vector_dofs);
    fluid_dofs.GetSubVector(vector_dofs, xval);
    state_in_current_element.UseExternalData(xval.GetData(), num_dof, num_equations);
    flux_in_current_element.SetSize(num_equations, dim*num_dof);
    dg_euler_operator.computeFluxVectorInElement(element, state_in_current_element, flux_in_current_element);

    for (int equation = 0; equation<num_equations; equation++) {
      flux_in_current_element.GetRow(equation, tempval);
      flux_for_current_equation.UseExternalData(tempval.GetData(),dim,num_dof);
      for (int idim=0; idim<dim; idim++) {
        double deltaij = (idim == equation - 1) ? 1.0 : 0.0;
        for (int jdof=0; jdof<num_dof; jdof++) {
          if (equation == 0) {
            EXPECT_EQ(flux_for_current_equation(idim,jdof),rho*velocity_field(idim));
          }
          else if (equation < 4) {
            EXPECT_EQ(flux_for_current_equation(idim,jdof),rho*velocity_field(equation-1)*velocity_field(idim) + deltaij*p);
          }
          else {
            EXPECT_EQ(flux_for_current_equation(idim,jdof),enthalpy*rho*velocity_field(idim));
          }
        }
      }
    }
  }
}

TEST(DGEulerAssembly, EulerFluxIsCorrectForSpatiallyVaryingStateIn3D) {
  int order = 1;
  int dim = 3;
  int num_equations = 5;
  double specific_heat_ratio = 1.4;
  const Species species{.specific_heat_ratio = specific_heat_ratio};
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(3, 3, 3, mfem::Element::HEXAHEDRON, 1, 1, 1);
  Discretization discretization(&mesh, order, FETypes::DG, num_equations);
  mfem::FiniteElementSpace finite_element_space = discretization.getFeSpace();
  DGEulerAssembly dg_euler_operator(finite_element_space, species);

  double tol = 1e-9;
  constexpr mfem::real_t c0(12.7), c1(-9.4), c2(2.2), c3(9.1);
  auto solution_vec = [&](const mfem::Vector &x, mfem::Vector &y) {
    mfem::real_t base_val = c0 + c1 * x[0] + c2 * x[1] + c3 * x[2];
    double rho = 0.5*x[0] + 1.0;
    double u = 2.0 * base_val;
    double v = 3.0 * base_val;
    double w = 4.0 * base_val;
    double p = 4.0 * base_val;
    double e = p/((specific_heat_ratio - 1.0)*rho)  + 0.5 * (u*u + v*v + w*w);
    y[0] = rho;
    y[1] = rho * u;
    y[2] = rho * v;
    y[3] = rho * w;
    y[4] = rho * e;
  };

  mfem::VectorFunctionCoefficient fluid_coeff(num_equations,solution_vec);
  mfem::GridFunction fluid_dofs(&discretization.getFeSpace());
  fluid_dofs.ProjectCoefficient(fluid_coeff);

  mfem::DenseMatrix flux_in_current_element;
  mfem::DenseMatrix state_in_current_element;
  mfem::Vector state_at_current_dof(num_equations);
  mfem::Array<int> vector_dofs;
  mfem::Vector xval,tempval;
  mfem::DenseMatrix flux_for_current_equation;

  for (int element=0; element<finite_element_space.GetNE(); element++) {
    int num_dof = finite_element_space.GetFE(element)->GetDof();
    finite_element_space.GetElementVDofs(element, vector_dofs);
    fluid_dofs.GetSubVector(vector_dofs, xval);
    state_in_current_element.UseExternalData(xval.GetData(), num_dof, num_equations);
    flux_in_current_element.SetSize(num_equations, dim*num_dof);
    dg_euler_operator.computeFluxVectorInElement(element, state_in_current_element, flux_in_current_element);

    const mfem::FiniteElement *fe = finite_element_space.GetFE(element);
    mfem::ElementTransformation *trans = finite_element_space.GetElementTransformation(element);
    const mfem::IntegrationRule &nodes = fe->GetNodes();

    for (int equation = 0; equation<num_equations; equation++) {
      flux_in_current_element.GetRow(equation, tempval);
      flux_for_current_equation.UseExternalData(tempval.GetData(),dim,num_dof);
      for (int idim=0; idim<dim; idim++) {
        double deltaij = (idim == equation - 1) ? 1.0 : 0.0;
        for (int jdof=0; jdof<num_dof; jdof++) {
          mfem::IntegrationPoint local_coords = nodes.IntPoint(jdof);
          mfem::Vector real_coords(dim);
          trans->Transform(local_coords, real_coords);
          mfem::Vector exact_solution_at_dof(5);
          solution_vec(real_coords,exact_solution_at_dof);
          double rho = exact_solution_at_dof[0];
          mfem::Vector velocity_field{exact_solution_at_dof[1]/rho,exact_solution_at_dof[2]/rho,exact_solution_at_dof[3]/rho};
          double p = (exact_solution_at_dof[4]/rho - 0.5 * (velocity_field[0]*velocity_field[0] + velocity_field[1]*velocity_field[1] + velocity_field[2]*velocity_field[2]))*((specific_heat_ratio - 1.0)*rho);
          double enthalpy = exact_solution_at_dof[4]/rho + p/rho;

          if (equation == 0) {
            EXPECT_NEAR(flux_for_current_equation(idim,jdof),rho*velocity_field(idim),tol);
          }
          else if (equation < 4) {
            EXPECT_NEAR(flux_for_current_equation(idim,jdof),rho*velocity_field(equation-1)*velocity_field(idim) + deltaij*p,tol);
          }
          else {
            EXPECT_NEAR(flux_for_current_equation(idim,jdof),enthalpy*rho*velocity_field(idim),tol);
          }
        }
      }
    }
  }
}

TEST(DGEulerOperator, EulerWeakDivergenceIsZeroForConstantStateIn3D) {
  int order = 1;
  int num_equations = 5;
  int dim = 3;
  double specific_heat_ratio = 1.4;
  const Species species{.specific_heat_ratio = specific_heat_ratio};
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(3, 3, 3, mfem::Element::HEXAHEDRON, 1, 1, 1);
  Discretization discretization(&mesh, order, FETypes::DG, num_equations);
  mfem::FiniteElementSpace finite_element_space = discretization.getFeSpace();
  DGEulerAssembly dg_euler_operator(finite_element_space, species);

  double rho = 1.0;
  double u = 2.0;
  double v = 3.0;
  double w = 4.0;
  double p = 5.0;
  mfem::Vector velocity_field{u,v,w};
  double e = p/((specific_heat_ratio - 1.0)*rho)  + 0.5 * (u*u + v*v + w*w);

  const mfem::Vector fluid_constant_vals{rho,rho*u,rho*v,rho*w,rho*e};
  mfem::VectorConstantCoefficient fluid_constant_coeff(fluid_constant_vals);
  mfem::GridFunction fluid_dofs(&finite_element_space);
  fluid_dofs.ProjectCoefficient(fluid_constant_coeff);

  mfem::Vector rhs(discretization.getFeSpace().GetTrueVSize());
  EXPECT_TRUE(rhs.Size() == fluid_dofs.Size());
  rhs = 0.0;

  mfem::Array<int> vector_dofs;
  mfem::DenseMatrix flux_in_current_element;
  mfem::DenseMatrix state_in_current_element;
  mfem::Vector vector_state_in_current_element;
  mfem::DenseMatrix rhs_in_current_element;
  mfem::Vector vector_rhs_in_current_element;

  for (int element=0; element<finite_element_space.GetNE(); element++)
  {
    int num_dof = finite_element_space.GetFE(element)->GetDof();
    finite_element_space.GetElementVDofs(element, vector_dofs);

    fluid_dofs.GetSubVector(vector_dofs, vector_state_in_current_element);
    state_in_current_element.UseExternalData(vector_state_in_current_element.GetData(), num_dof, num_equations);
    flux_in_current_element.SetSize(num_equations, dim*num_dof);
    dg_euler_operator.computeFluxVectorInElement(element,state_in_current_element,flux_in_current_element);

    rhs.GetSubVector(vector_dofs, vector_rhs_in_current_element);
    rhs_in_current_element.UseExternalData(vector_rhs_in_current_element.GetData(), num_dof, num_equations);
    flux_in_current_element.SetSize(num_equations, dim*num_dof);
    dg_euler_operator.applyWeakDivergenceInElement(element, flux_in_current_element, rhs_in_current_element);

    rhs.SetSubVector(vector_dofs, vector_rhs_in_current_element);

  }

  auto rhs_l2  = rhs.Sum();
  EXPECT_NEAR(rhs_l2, 0., 1e-12);
}

TEST(DGEulerOperator, EulerWeakDivergenceIsZeroForOrderZeroIn3D) {
  int order = 0;
  int num_equations = 5;
  int dim = 3;
  double specific_heat_ratio = 1.4;
  const Species species{.specific_heat_ratio = specific_heat_ratio};
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(3, 3, 3, mfem::Element::HEXAHEDRON, 1, 1, 1);
  Discretization discretization(&mesh, order, FETypes::DG, num_equations);
  mfem::FiniteElementSpace finite_element_space = discretization.getFeSpace();
  DGEulerAssembly dg_euler_operator(finite_element_space, species);

  auto solution_vec = [&](const mfem::Vector &x, mfem::Vector &y) {
    double rho = sin(x[0])+1.0;
    double u = cos(x[0])+2.0;
    double v = sin(x[1])+3.0;
    double w = cos(x[1])+4.0;
    double p = sin(x[2])+5.0;
    double e = p/((specific_heat_ratio - 1.0)*rho)  + 0.5 * (u*u + v*v + w*w);
    y[0] = rho;
    y[1] = rho * u;
    y[2] = rho * v;
    y[3] = rho * w;
    y[4] = rho * e;
  };

  mfem::VectorFunctionCoefficient fluid_coeff(num_equations,solution_vec);
  mfem::GridFunction fluid_dofs(&discretization.getFeSpace());
  fluid_dofs.ProjectCoefficient(fluid_coeff);

  mfem::Vector rhs(discretization.getFeSpace().GetTrueVSize());
  EXPECT_TRUE(rhs.Size() == fluid_dofs.Size());
  rhs = 0.0;
  mfem::Array<int> vector_dofs;
  mfem::DenseMatrix flux_in_current_element;
  mfem::DenseMatrix state_in_current_element;
  mfem::Vector vector_state_in_current_element;
  mfem::DenseMatrix rhs_in_current_element;
  mfem::Vector vector_rhs_in_current_element;

  for (int element=0; element<finite_element_space.GetNE(); element++)
  {
    int num_dof = finite_element_space.GetFE(element)->GetDof();
    finite_element_space.GetElementVDofs(element, vector_dofs);

    fluid_dofs.GetSubVector(vector_dofs, vector_state_in_current_element);
    state_in_current_element.UseExternalData(vector_state_in_current_element.GetData(), num_dof, num_equations);
    flux_in_current_element.SetSize(num_equations, dim*num_dof);
    dg_euler_operator.computeFluxVectorInElement(element,state_in_current_element,flux_in_current_element);

    rhs.GetSubVector(vector_dofs, vector_rhs_in_current_element);
    rhs_in_current_element.UseExternalData(vector_rhs_in_current_element.GetData(), num_dof, num_equations);
    flux_in_current_element.SetSize(num_equations, dim*num_dof);
    dg_euler_operator.applyWeakDivergenceInElement(element, flux_in_current_element, rhs_in_current_element);

    rhs.SetSubVector(vector_dofs, vector_rhs_in_current_element);

  }

  auto rhs_l2  = rhs.Norml2();
  EXPECT_NEAR(rhs_l2, 0., 1e-13);
}

TEST(DGAssembly, ApplyInverseMassMatrixRecoversProjectionCoefficients) {
  int order = 1;
  int num_equations = 5;
  double specific_heat_ratio = 1.4;
  const Species species{.specific_heat_ratio = specific_heat_ratio};
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(3, 3, 3, mfem::Element::HEXAHEDRON, 1, 1, 1);
  Discretization discretization(&mesh, order, FETypes::DG, num_equations);
  mfem::FiniteElementSpace finite_element_space = discretization.getFeSpace();
  DGEulerAssembly dg_euler_operator(finite_element_space, species);

  auto solution_vec = [&](const mfem::Vector &x, mfem::Vector &y) {
    double rho = sin(x[0])+1.0;
    double u = cos(x[0])+2.0;
    double v = sin(x[1])+3.0;
    double w = cos(x[1])+4.0;
    double p = sin(x[2])+5.0;
    double e = p/((specific_heat_ratio - 1.0)*rho)  + 0.5 * (u*u + v*v + w*w);
    y[0] = rho;
    y[1] = rho * u;
    y[2] = rho * v;
    y[3] = rho * w;
    y[4] = rho * e;
  };

  mfem::VectorFunctionCoefficient fluid_coeff(num_equations,solution_vec);

  mfem::GridFunction c_true(&discretization.getFeSpace());
  c_true.ProjectCoefficient(fluid_coeff);

  mfem::LinearForm b(&finite_element_space);
  b.AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(fluid_coeff));
  b.Assemble();

  mfem::Vector c_test(finite_element_space.GetVSize());
  dg_euler_operator.applyInverseMass(b,c_test);

  mfem::Vector diff = c_true;
  auto exact_l2 = diff.Norml2();
  diff.Add(-1.0, c_test);
  auto diff_l2 = diff.Norml2();

  EXPECT_NEAR(diff_l2/exact_l2, 0., 1e-13);

};

TEST(DGEulerOperator, EulerBoundaryFluxIsZeroForConstantStateIn3D) {
  int order = 1;
  int num_equations = 5;
  double specific_heat_ratio = 1.4;
  const Species species{.specific_heat_ratio = specific_heat_ratio};
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(3, 3, 3, mfem::Element::HEXAHEDRON, 1, 1, 1);
  Discretization discretization(&mesh, order, FETypes::DG, num_equations);
  mfem::FiniteElementSpace finite_element_space = discretization.getFeSpace();
  DGEulerAssembly dg_euler_operator(finite_element_space, species);

  double rho = 2.0;
  double u = 3.0;
  double v = 4.0;
  double w = 5.0;
  double p = 6.0;
  double e = p/((specific_heat_ratio - 1.0)*rho)  + 0.5 * (u*u + v*v + w*w);
  const mfem::Vector fluid_constant_vals{rho,rho*u,rho*v,rho*w,rho*e};
  mfem::VectorConstantCoefficient fluid_constant_coeff(fluid_constant_vals);
  mfem::GridFunction fluid_dofs(&finite_element_space);
  fluid_dofs.ProjectCoefficient(fluid_constant_coeff);

  mfem::Vector rhs(discretization.getFeSpace().GetTrueVSize());
  EXPECT_TRUE(rhs.Size() == fluid_dofs.Size());
  rhs = 0.0;
  dg_euler_operator.computeFluxOnElementBoundaries(fluid_dofs,rhs);

  auto rhs_sum  = rhs.Sum();
  EXPECT_NEAR(rhs_sum, 0., 1e-12);
}

TEST(DGEulerAssembly, checkNoWork) {

  // check that the dg basis cross rt dot dg = 0
  // this implies v cross b dot v = 0, or that there is no work done by
  // that portion of the lorentz force

  constexpr double absolute_tolerance = 1e-13;
  // single element, there's probably a better way to do this
  constexpr int nx = 1;
  constexpr mfem::real_t length = 1.0;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(nx, nx, nx, mfem::Element::HEXAHEDRON, length, length, length);

  constexpr int basis_order = 0;
  constexpr int dim = 3;

  mfem::DG_FECollection dg_finite_element_collection(basis_order, dim);
  mfem::FiniteElementSpace dg_finite_element_space(&mesh, &dg_finite_element_collection, 1, mfem::Ordering::byNODES);

  mfem::RT_FECollection rt_finite_element_collection(basis_order + 1, dim);
  mfem::FiniteElementSpace rt_finite_element_space(&mesh, &rt_finite_element_collection, 1, mfem::Ordering::byNODES);

  const int ielem = 0;

  mfem::Vector t1{.82,.04,0.7};
  mfem::Vector t2{.126,.88,0.13};
  mfem::Vector t3{.53,.65,0.23};
  mfem::Vector t4{.99,.11,0.45};

  std::vector<mfem::Vector> test_points{t1,t2,t3,t4};

  for (const auto & point : test_points) {

    auto elem = mesh.GetElement(ielem);

    mfem::Array<int> vertices;
    elem->GetVertices(vertices);

    mfem::IsoparametricTransformation element_transformation;
    mesh.GetElementTransformation(ielem, &element_transformation);

    mfem::IntegrationPoint integration_point;
    element_transformation.TransformBack(point, integration_point);
    element_transformation.SetIntPoint(&integration_point);

    auto dg_finite_element = dg_finite_element_space.GetFE(ielem);
    auto rt_finite_element = rt_finite_element_space.GetFE(ielem);
    mfem::Vector dg_basis_vals(dg_finite_element->GetDof());
    mfem::DenseMatrix rt_basis_vals(rt_finite_element->GetDof(),dim);
    dg_finite_element->CalcShape(integration_point, dg_basis_vals);
    rt_finite_element->CalcVShape(integration_point, rt_basis_vals);

    mfem::Vector dg_at_point{0.,0.,0.};
    mfem::Vector rt_at_point{0.,0.,0.};

    for (int idof = 0; idof < dg_finite_element->GetDof(); ++idof) {
      for (int idim = 0; idim < 3; ++idim) {
        dg_at_point(idim) += element_transformation.Weight() * dg_basis_vals(idof);
      }
    }

    for (int idof = 0; idof < rt_finite_element->GetDof(); ++idof) {
      for (int idim = 0; idim < 3; ++idim) {
        rt_at_point(idim) += element_transformation.Weight() * rt_basis_vals(idof,idim);
      }
    }

    mfem::Vector cross{0.,0.,0.};
    dg_at_point.cross3D(rt_at_point, cross);
    const double cross_dot = cross * dg_at_point;

    EXPECT_NEAR(cross_dot, 0., absolute_tolerance);

  }

}

TEST(DGEulerAssembly, electrostaticSourceTermOrder0) {
  constexpr double tolerance = 1e-13;
  constexpr int nx = 5;
  constexpr mfem::real_t length = 1.0;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(nx, nx, nx, mfem::Element::HEXAHEDRON, length, length, length);

  constexpr int basis_order = 0;
  constexpr int num_equations = 5;

  constexpr double charge_over_mass = 8.92;

  Discretization fluid_discretization(&mesh, basis_order, FETypes::DG, num_equations);
  Discretization potential_discretization(&mesh, basis_order + 1, FETypes::HGRAD);

  const mfem::Vector fluid_constant_vals{1.2,-3.5,4.7,2.1,12.8};
  mfem::VectorConstantCoefficient fluid_constant_coeff(fluid_constant_vals);

  mfem::GridFunction fluid_state(&fluid_discretization.getFeSpace());
  fluid_state.ProjectCoefficient(fluid_constant_coeff);

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

  for (const auto & add_energy_source : {true, false})
  {
    mfem::LinearForm lin_form(&fluid_discretization.getFeSpace());
    lin_form.AddDomainIntegrator(new EulerMaxwellSourceIntegrator(fluid_state, es_field_state, charge_over_mass, add_energy_source));
    lin_form.Assemble();

    mfem::Vector source_exact{
      0.,
      charge_over_mass * fluid_constant_vals[0] * e[0],
      charge_over_mass * fluid_constant_vals[0] * e[1],
      charge_over_mass * fluid_constant_vals[0] * e[2],
      add_energy_source * charge_over_mass * (
        fluid_constant_vals[1] * e[0] + fluid_constant_vals[2] * e[1] + fluid_constant_vals[3] * e[2]
      )
    };
    mfem::VectorConstantCoefficient source_exact_coeff(source_exact);
    mfem::LinearForm exact_form(&fluid_discretization.getFeSpace());
    exact_form.AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(source_exact_coeff));
    exact_form.Assemble();

    mfem::Vector diff = exact_form;
    auto exact_l2 = diff.Norml2();
    diff.Add(-1.0, lin_form);
    auto diff_l2 = diff.Norml2();

    EXPECT_NEAR(diff_l2/exact_l2, 0., tolerance);
  }

}

TEST(DGEulerAssembly, electrostaticSourceTermOrder1) {
  constexpr double tolerance = 1e-13;
  constexpr int nx = 5;
  constexpr mfem::real_t length = 1.0;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(nx, nx, nx, mfem::Element::HEXAHEDRON, length, length, length);

  constexpr int basis_order = 1;
  constexpr int num_equations = 5;

  constexpr double charge_over_mass = 8.92;

  Discretization fluid_discretization(&mesh, basis_order, FETypes::DG, num_equations);
  Discretization potential_discretization(&mesh, basis_order, FETypes::HGRAD);

  constexpr mfem::real_t c0(12.7), c1(-9.4), c2(2.2), c3(9.1);
  auto linear_vec = [&](const mfem::Vector &x, mfem::Vector &y) {
    mfem::real_t base_val = c0 + c1 * x[0] + c2 * x[1] + c3 * x[2];
    for (int i = 0; i < num_equations; ++i) {
      y[i] = (i+1) * base_val;
    }
  };
  mfem::VectorFunctionCoefficient fluid_coeff(num_equations,linear_vec);

  mfem::GridFunction fluid_state(&fluid_discretization.getFeSpace());
  fluid_state.ProjectCoefficient(fluid_coeff);

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

  for (const auto & add_energy_source : {true, false})
  {
    mfem::LinearForm lin_form(&fluid_discretization.getFeSpace());
    lin_form.AddDomainIntegrator(new EulerMaxwellSourceIntegrator(fluid_state, es_field_state, charge_over_mass, add_energy_source));
    lin_form.Assemble();

    auto source_exact_vec = [&](const mfem::Vector &x, mfem::Vector &y) {
      linear_vec(x,y);
      const auto density = y[0];
      const mfem::Vector velocity{y[1]/density, y[2]/density, y[3]/density};
      y[0] = 0;
      y[1] = charge_over_mass * density * e(0);
      y[2] = charge_over_mass * density * e(1);
      y[3] = charge_over_mass * density * e(2);
      y[4] = add_energy_source ? (charge_over_mass * density * (e * velocity)) : 0.;
    };
    mfem::VectorFunctionCoefficient source_exact_coeff(num_equations,source_exact_vec);
    mfem::LinearForm exact_form(&fluid_discretization.getFeSpace());
    exact_form.AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(source_exact_coeff));
    exact_form.Assemble();

    mfem::Vector diff = exact_form;
    auto exact_l2 = diff.Norml2();
    diff.Add(-1.0, lin_form);
    auto diff_l2 = diff.Norml2();

    EXPECT_NEAR(diff_l2/exact_l2, 0., tolerance);
  }

}

TEST(DGEulerAssembly, computeSourcesElectrostatic) {
  constexpr double tolerance = 1e-13;
  constexpr int nx = 5;
  constexpr mfem::real_t length = 1.0;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(nx, nx, nx, mfem::Element::HEXAHEDRON, length, length, length);

  constexpr int basis_order = 1;
  constexpr int num_equations = 5;

  constexpr double charge_over_mass = 8.92;
  const Species species{.charge_over_mass = charge_over_mass};

  Discretization fluid_discretization(&mesh, basis_order, FETypes::DG, num_equations);
  Discretization potential_discretization(&mesh, basis_order, FETypes::HGRAD);

  constexpr mfem::real_t c0(12.7), c1(-9.4), c2(2.2), c3(9.1);
  auto linear_vec = [&](const mfem::Vector &x, mfem::Vector &y) {
    mfem::real_t base_val = c0 + c1 * x[0] + c2 * x[1] + c3 * x[2];
    for (int i = 0; i < num_equations; ++i) {
      y[i] = (i+1) * base_val;
    }
  };
  mfem::VectorFunctionCoefficient fluid_coeff(num_equations,linear_vec);

  LowFidelitySpeciesState fluid_state(fluid_discretization, species);
  mfem::GridFunction& fluid_grid_function = fluid_state.getGridFunction();
  fluid_grid_function.ProjectCoefficient(fluid_coeff);

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

  DGEulerAssembly op(fluid_discretization.getFeSpace(), species);

  mfem::Vector rhs(fluid_discretization.getFeSpace().GetTrueVSize());
  constexpr double rhs_offset = 123.456;

  rhs += rhs_offset;

  EXPECT_TRUE(rhs.Size() == fluid_grid_function.Size());
  op.computeSources(fluid_state, es_field_state, rhs);

  auto source_exact_vec = [&](const mfem::Vector &x, mfem::Vector &y) {
    linear_vec(x,y);
    const auto density = y[0];
    const mfem::Vector velocity{y[1]/density, y[2]/density, y[3]/density};
    y[0] = 0;
    y[1] = charge_over_mass * density * e(0);
    y[2] = charge_over_mass * density * e(1);
    y[3] = charge_over_mass * density * e(2);
    y[4] = 0.;
  };
  mfem::VectorFunctionCoefficient source_exact_coeff(num_equations,source_exact_vec);
  mfem::LinearForm exact_form(&fluid_discretization.getFeSpace());
  exact_form.AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(source_exact_coeff));
  exact_form.Assemble();

  mfem::Vector diff = exact_form;
  diff += rhs_offset;
  auto exact_l2 = diff.Norml2();
  diff.Add(-1.0, rhs);
  auto diff_l2 = diff.Norml2();

  EXPECT_NEAR(diff_l2/exact_l2, 0., tolerance);
}

TEST(DGEulerAssembly, KineticEnergyIntegratorOrder1) {
  constexpr double tolerance = 1e-13;
  constexpr int nx = 5;
  constexpr mfem::real_t length = 1.0;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(nx, nx, nx, mfem::Element::HEXAHEDRON, length, length, length);

  constexpr int basis_order = 1;
  constexpr int num_equations = 5;

  Discretization fluid_discretization(&mesh, basis_order, FETypes::DG, num_equations);

  constexpr mfem::real_t c0(12.7), c1(-9.4), c2(2.2), c3(9.1);
  auto linear_vec = [&](const mfem::Vector &x, mfem::Vector &y) {
    mfem::real_t base_val = c0 + c1 * x[0] + c2 * x[1] + c3 * x[2];
    for (int i = 0; i < num_equations; ++i) {
      y[i] = (i+1) * base_val;
    }
  };
  mfem::VectorFunctionCoefficient fluid_coeff(num_equations,linear_vec);

  mfem::GridFunction fluid_state(&fluid_discretization.getFeSpace());
  fluid_state.ProjectCoefficient(fluid_coeff);

  mfem::LinearForm lin_form(&fluid_discretization.getFeSpace());
  lin_form.AddDomainIntegrator(new EulerKineticEnergyIntegrator(fluid_state));
  lin_form.Assemble();

  auto source_exact_vec = [&](const mfem::Vector &x, mfem::Vector &y) {
    linear_vec(x,y);
    const auto density = y[0];
    const mfem::Vector velocity{y[1]/density, y[2]/density, y[3]/density};
    const mfem::Vector momentum{y[1], y[2], y[3]};
    y[0] = 0;
    y[1] = 0;
    y[2] = 0;
    y[3] = 0;
    y[4] = 1./2. * (velocity * momentum);
  };
  mfem::VectorFunctionCoefficient source_exact_coeff(num_equations,source_exact_vec);
  mfem::LinearForm exact_form(&fluid_discretization.getFeSpace());
  exact_form.AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(source_exact_coeff));
  exact_form.Assemble();

  mfem::Vector diff = exact_form;
  auto exact_l2 = diff.Norml2();
  diff.Add(-1.0, lin_form);
  auto diff_l2 = diff.Norml2();

  EXPECT_NEAR(diff_l2/exact_l2, 0., tolerance);

}

TEST(DGEulerAssembly, computeIntegratedKineticEnergy) {
  constexpr double tolerance = 1e-13;
  constexpr int nx = 5;
  constexpr mfem::real_t length = 1.0;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(nx, nx, nx, mfem::Element::HEXAHEDRON, length, length, length);

  constexpr int basis_order = 1;
  constexpr int num_equations = 5;

  constexpr double charge_over_mass = 8.92;
  const Species species{.charge_over_mass = charge_over_mass};

  Discretization fluid_discretization(&mesh, basis_order, FETypes::DG, num_equations);

  constexpr mfem::real_t c0(12.7), c1(-9.4), c2(2.2), c3(9.1);
  auto linear_vec = [&](const mfem::Vector &x, mfem::Vector &y) {
    mfem::real_t base_val = c0 + c1 * x[0] + c2 * x[1] + c3 * x[2];
    for (int i = 0; i < num_equations; ++i) {
      y[i] = (i+1) * base_val;
    }
  };
  mfem::VectorFunctionCoefficient fluid_coeff(num_equations,linear_vec);

  LowFidelitySpeciesState fluid_state(fluid_discretization, species);
  mfem::GridFunction& fluid_grid_function = fluid_state.getGridFunction();
  fluid_grid_function.ProjectCoefficient(fluid_coeff);

  DGEulerAssembly op(fluid_discretization.getFeSpace(), species);

  mfem::Vector rhs(fluid_discretization.getFeSpace().GetTrueVSize());
  constexpr double rhs_offset = 123.456;

  rhs += rhs_offset;

  EXPECT_TRUE(rhs.Size() == fluid_grid_function.Size());
  op.computeIntegratedKineticEnergy(fluid_state, rhs);

  auto source_exact_vec = [&](const mfem::Vector &x, mfem::Vector &y) {
    linear_vec(x,y);
    const auto density = y[0];
    const mfem::Vector velocity{y[1]/density, y[2]/density, y[3]/density};
    const mfem::Vector momentum{y[1], y[2], y[3]};
    y[0] = 0;
    y[1] = 0;
    y[2] = 0;
    y[3] = 0;
    y[4] = 1./2. * (velocity * momentum);
  };
  mfem::VectorFunctionCoefficient source_exact_coeff(num_equations,source_exact_vec);
  mfem::LinearForm exact_form(&fluid_discretization.getFeSpace());
  exact_form.AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(source_exact_coeff));
  exact_form.Assemble();

  mfem::Vector diff = exact_form;
  diff += rhs_offset;
  auto exact_l2 = diff.Norml2();
  diff.Add(-1.0, rhs);
  auto diff_l2 = diff.Norml2();

  EXPECT_NEAR(diff_l2/exact_l2, 0., tolerance);
}

} // namespace
