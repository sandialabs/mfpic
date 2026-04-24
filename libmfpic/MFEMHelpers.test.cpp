#include <libmfpic/MFEMHelpers.hpp>

#include <libmfpic/Discretization.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

void checkTransformedVectorCoefficientEvaluatesToExpectedConstant(
  TransformedVectorCoefficient& transformed_vector_coefficient,
  mfem::Mesh& mesh,
  const double expected_value)
{
  Discretization identity_test_function_discretization = getIdentityTestFunctionDiscretization(mesh);

  mfem::LinearForm integrals_by_cell(&identity_test_function_discretization.getFeSpace());
  integrals_by_cell.AddDomainIntegrator(new mfem::DomainLFIntegrator(transformed_vector_coefficient));
  integrals_by_cell.Assemble();

  const double cell_volume = mesh.GetElementVolume(0);
  const double expected_integral = expected_value * cell_volume;
  constexpr double absolute_tolerance = 1e-13;
  for (int i = 0; i < integrals_by_cell.Size(); ++i){
    EXPECT_NEAR(expected_integral, integrals_by_cell[i], absolute_tolerance);
  }
}

TEST(MFEMHelpers, TransformedVectorCoefficient1D) {
  const mfem::Vector vec{1.2, 2.3, 3.4};
  auto vector_coefficient = std::make_unique<mfem::VectorConstantCoefficient>(vec);

  auto transformation = [](const mfem::Vector& vec) {
    return 4.5 * (vec * vec);
  };

  TransformedVectorCoefficient transformed_vector_coefficient(std::move(vector_coefficient), transformation);

  constexpr int num_elems = 10;
  constexpr double length = 2.0;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems, length);

  const double expected_value = transformation(vec);
  checkTransformedVectorCoefficientEvaluatesToExpectedConstant(transformed_vector_coefficient, mesh, expected_value);
}

TEST(MFEMHelpers, TransformedVectorCoefficient2D) {
  constexpr int vector_dim = 2;
  const mfem::Vector vec{5.6, 7.8};
  auto vector_function = [&](const mfem::Vector&, mfem::Vector& y) { y = vec; };
  auto vector_coefficient = std::make_unique<mfem::VectorFunctionCoefficient>(vector_dim, vector_function);

  auto transformation = [](const mfem::Vector& vec) {
    return 2.3 * vec[0] + vec[1] * vec[1];
  };

  TransformedVectorCoefficient transformed_vector_coefficient(std::move(vector_coefficient), transformation);

  constexpr int num_elems_per_dim = 5;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(num_elems_per_dim, num_elems_per_dim, mfem::Element::QUADRILATERAL);

  const double expected_value = transformation(vec);
  checkTransformedVectorCoefficientEvaluatesToExpectedConstant(transformed_vector_coefficient, mesh, expected_value);
}

TEST(MFEMHelpers, TransformedVectorCoefficient3D) {
  const mfem::Vector vec{6.3, 1.6, 9.1};
  auto vector_coefficient = std::make_unique<mfem::VectorConstantCoefficient>(vec);

  auto transformation = [](const mfem::Vector& vec) {
    return sin(vec[0]) + cos(vec[1]) + exp(vec[2]);
  };

  TransformedVectorCoefficient transformed_vector_coefficient(std::move(vector_coefficient), transformation);

  constexpr int num_elems_per_dim = 7;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
    num_elems_per_dim, num_elems_per_dim, num_elems_per_dim, mfem::Element::TETRAHEDRON);

  const double expected_value = transformation(vec);
  checkTransformedVectorCoefficientEvaluatesToExpectedConstant(transformed_vector_coefficient, mesh, expected_value);
}

}