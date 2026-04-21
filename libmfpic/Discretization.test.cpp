#include <libmfpic/Discretization.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

void checkIdentityTestFunctionDiscretizationIntegratesConstant(Discretization& identity_test_function_discretization) {
  constexpr double constant = 3.14;
  mfem::ConstantCoefficient constant_coefficient(constant);
  mfem::FiniteElementSpace* finite_element_space = &identity_test_function_discretization.getFeSpace();
  mfem::LinearForm elementwise_integrals(finite_element_space);
  elementwise_integrals.AddDomainIntegrator(new mfem::DomainLFIntegrator(constant_coefficient));

  elementwise_integrals.Assemble();

  mfem::Mesh* mesh = finite_element_space->GetMesh();
  const double delta_x = mesh->GetElementSize(0);
  const double expected_element_integral = constant * delta_x;
  for (int i = 0; i < elementwise_integrals.Size(); ++i) {
    EXPECT_DOUBLE_EQ(expected_element_integral, elementwise_integrals[i]);
  }
}

TEST(Discretization, IdentityTestFunctionDiscretizationFromMesh) {
  constexpr int num_elems = 10;
  constexpr double length = 2.0;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems, length);

  Discretization identity_test_function_discretization = getIdentityTestFunctionDiscretization(mesh);
  checkIdentityTestFunctionDiscretizationIntegratesConstant(identity_test_function_discretization);
}

TEST(Discretization, IdentityTestFunctionDiscretizationFromDiscretization) {
  constexpr int num_elems = 10;
  constexpr double length = 2.0;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems, length);

  constexpr int basis_order = 1;
  Discretization discretization(&mesh, basis_order, FETypes::HGRAD);

  Discretization identity_test_function_discretization = getIdentityTestFunctionDiscretization(discretization);
  checkIdentityTestFunctionDiscretizationIntegratesConstant(identity_test_function_discretization);
}

}