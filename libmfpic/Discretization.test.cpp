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
  const double cell_volume = mesh->GetElementVolume(0);
  const double expected_element_integral = constant * cell_volume;
  constexpr double absolute_tolerance = 1e-15;
  for (int i = 0; i < elementwise_integrals.Size(); ++i) {
    EXPECT_NEAR(expected_element_integral, elementwise_integrals[i], absolute_tolerance);
  }
}

TEST(Discretization, IdentityTestFunctionDiscretizationFromMesh1D) {
  constexpr int num_elems = 10;
  constexpr double length = 2.0;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems, length);

  Discretization identity_test_function_discretization = getIdentityTestFunctionDiscretization(mesh);
  checkIdentityTestFunctionDiscretizationIntegratesConstant(identity_test_function_discretization);
}

TEST(Discretization, IdentityTestFunctionDiscretizationFromMesh2D) {
  constexpr int num_elems_per_dim = 11;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(num_elems_per_dim, num_elems_per_dim, mfem::Element::QUADRILATERAL);

  Discretization identity_test_function_discretization = getIdentityTestFunctionDiscretization(mesh);
  checkIdentityTestFunctionDiscretizationIntegratesConstant(identity_test_function_discretization);
}

TEST(Discretization, IdentityTestFunctionDiscretizationFromMesh3D) {
  constexpr int num_elems_per_dim = 12;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
    num_elems_per_dim, num_elems_per_dim, num_elems_per_dim, mfem::Element::HEXAHEDRON);

  Discretization identity_test_function_discretization = getIdentityTestFunctionDiscretization(mesh);
  checkIdentityTestFunctionDiscretizationIntegratesConstant(identity_test_function_discretization);
}

TEST(Discretization, IdentityTestFunctionDiscretizationFromDiscretization1D) {
  constexpr int num_elems = 11;
  constexpr double length = 1.5;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems, length);

  constexpr int basis_order = 1;
  Discretization discretization(&mesh, basis_order, FETypes::HGRAD);

  Discretization identity_test_function_discretization = getIdentityTestFunctionDiscretization(discretization);
  checkIdentityTestFunctionDiscretizationIntegratesConstant(identity_test_function_discretization);
}

TEST(Discretization, IdentityTestFunctionDiscretizationFromDiscretization2D) {
  constexpr int num_elems_per_dim = 9;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(num_elems_per_dim, num_elems_per_dim, mfem::Element::TRIANGLE);

  constexpr int basis_order = 0;
  constexpr int num_fields = 2;
  Discretization discretization(&mesh, basis_order, FETypes::DG, num_fields);

  Discretization identity_test_function_discretization = getIdentityTestFunctionDiscretization(discretization);
  checkIdentityTestFunctionDiscretizationIntegratesConstant(identity_test_function_discretization);
}

TEST(Discretization, IdentityTestFunctionDiscretizationFromDiscretization3D) {
  constexpr int num_elems_per_dim = 8;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
    num_elems_per_dim, num_elems_per_dim, num_elems_per_dim, mfem::Element::TETRAHEDRON);

  constexpr int basis_order = 2;
  Discretization discretization(&mesh, basis_order, FETypes::HGRAD);

  Discretization identity_test_function_discretization = getIdentityTestFunctionDiscretization(discretization);
  checkIdentityTestFunctionDiscretizationIntegratesConstant(identity_test_function_discretization);
}

}