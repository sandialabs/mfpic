#include <libmfpic/MFEMHelpers.hpp>

#include <libmfpic/Discretization.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

TEST(MFEMHelpers, TransformedVectorCoefficient) {
  const mfem::Vector vec{1.2, 2.3, 3.4};
  auto vector_coefficient = std::make_unique<mfem::VectorConstantCoefficient>(vec);

  auto transformation = [](const mfem::Vector& vec) {
    return 4.5 * (vec * vec);
  };

  TransformedVectorCoefficient transformed_vector_coefficient(std::move(vector_coefficient), transformation);

  constexpr int num_elems = 10;
  constexpr double length = 2.0;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems, length);
  Discretization identity_test_function_discretization = getIdentityTestFunctionDiscretization(mesh);

  mfem::LinearForm integrals_by_cell(&identity_test_function_discretization.getFeSpace());
  integrals_by_cell.AddDomainIntegrator(new mfem::DomainLFIntegrator(transformed_vector_coefficient));
  integrals_by_cell.Assemble();

  const double expected_integrand = transformation(vec);
  constexpr double dx = length / num_elems;
  const double expected_integral = expected_integrand * dx;
  for (int i = 0; i < integrals_by_cell.Size(); ++i){
    EXPECT_DOUBLE_EQ(expected_integral, integrals_by_cell[i]);
  }
}

}