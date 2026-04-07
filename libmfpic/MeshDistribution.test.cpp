#include <libmfpic/Errors.hpp>
#include <libmfpic/MeshDistribution.hpp>
#include <libmfpic/MeshUtilities.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

double uniformDistribution(const mfem::Vector&) {
  return 1.0;
}

void testThatUniformlyRandomlyGeneratedElementsAreActuallyDistributedUniformly(
  std::shared_ptr<mfem::Mesh> mesh,
  int num_points_to_generate = 20000,
  double relative_tolerance = 0.1
) {
  const int num_elements = mesh->GetNE();
  std::vector<int> num_points_generated_in_element(num_elements, 0);

  std::mt19937 generator;
  MeshDistribution distribution(mesh, uniformDistribution);
  for (int i = 0; i < num_points_to_generate; i++) {
    int element;
    std::tie(std::ignore, element) = distribution.generateRandomPointAndElement(generator);
    num_points_generated_in_element[element] += 1;
  }

  const double expected_points_generated_per_element = static_cast<double>(num_points_to_generate) / num_elements;
  const double absolute_tolerance = expected_points_generated_per_element * relative_tolerance;
  for (int num_points_generated : num_points_generated_in_element) {
    EXPECT_NEAR(num_points_generated, expected_points_generated_per_element, absolute_tolerance);
  }
}

void testThatUniformlyRandomlyGeneratedPointsAreActuallyDistributedUniformlyInASingleDimension(
  std::shared_ptr<mfem::Mesh> mesh,
  int dimension,
  int num_points_to_generate = 20000,
  int num_bins = 5,
  double relative_tolerance = 0.1
) {
  ASSERT_LT(dimension, mesh->Dimension());
  std::vector<int> bins(num_bins, 0);

  std::mt19937 generator;
  MeshDistribution distribution(mesh, uniformDistribution);
  for (int i = 0; i < num_points_to_generate; i++) {
    mfem::Vector point;
    std::tie(point, std::ignore) = distribution.generateRandomPointAndElement(generator);
    bins[point[dimension] * num_bins] += 1;
  }

  const double expected_points_generated_per_bin = static_cast<double>(num_points_to_generate) / num_bins;
  const double absolute_tolerance = expected_points_generated_per_bin * relative_tolerance;
  for (int num_points_generated : bins) {
    EXPECT_NEAR(num_points_generated, expected_points_generated_per_bin, absolute_tolerance);
  }
}

void testThatUniformlyRandomlyGeneratedPointsAreActuallyDistributedUniformly(
  std::shared_ptr<mfem::Mesh> mesh,
  int num_points_to_generate = 20000,
  int num_bins = 5,
  double relative_tolerance = 0.1
) {
  for (int dimension = 0; dimension < mesh->Dimension(); dimension++) {
    testThatUniformlyRandomlyGeneratedPointsAreActuallyDistributedUniformlyInASingleDimension(
      mesh,
      dimension,
      num_points_to_generate,
      num_bins,
      relative_tolerance
    );
  }
}

TEST(MeshDistribution, TestThatUniformlyRandomlyGeneratedPointsAreActuallyDistributedUniformlyFor1DMeshes) {
  constexpr int num_elements = 10;
  auto mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elements));

  testThatUniformlyRandomlyGeneratedElementsAreActuallyDistributedUniformly(mesh);
  testThatUniformlyRandomlyGeneratedPointsAreActuallyDistributedUniformly(mesh);
}

TEST(MeshDistribution, TestThatUniformlyRandomlyGeneratedPointsAreActuallyDistributedUniformlyForQuadMeshes) {
  auto mesh = std::make_shared<mfem::Mesh>(createMeshOfUnitBoxWith2ElemsPerDimension(mfem::Element::QUADRILATERAL));

  testThatUniformlyRandomlyGeneratedElementsAreActuallyDistributedUniformly(mesh);
  testThatUniformlyRandomlyGeneratedPointsAreActuallyDistributedUniformly(mesh);
}

TEST(MeshDistribution, TestThatUniformlyRandomlyGeneratedPointsAreActuallyDistributedUniformlyForTriMeshes) {
  auto mesh = std::make_shared<mfem::Mesh>(createMeshOfUnitBoxWith2ElemsPerDimension(mfem::Element::TRIANGLE));

  testThatUniformlyRandomlyGeneratedElementsAreActuallyDistributedUniformly(mesh);
  testThatUniformlyRandomlyGeneratedPointsAreActuallyDistributedUniformly(mesh);
}

TEST(MeshDistribution, TestThatUniformlyRandomlyGeneratedPointsAreActuallyDistributedUniformlyForHexMeshes) {
  auto mesh = std::make_shared<mfem::Mesh>(createMeshOfUnitBoxWith2ElemsPerDimension(mfem::Element::HEXAHEDRON));

  testThatUniformlyRandomlyGeneratedElementsAreActuallyDistributedUniformly(mesh);
  testThatUniformlyRandomlyGeneratedPointsAreActuallyDistributedUniformly(mesh);
}

TEST(MeshDistribution, TestThatUniformlyRandomlyGeneratedPointsAreActuallyDistributedUniformlyForTetMeshes) {
  auto mesh = std::make_shared<mfem::Mesh>(createMeshOfUnitBoxWith2ElemsPerDimension(mfem::Element::TETRAHEDRON));

  testThatUniformlyRandomlyGeneratedElementsAreActuallyDistributedUniformly(mesh);
  testThatUniformlyRandomlyGeneratedPointsAreActuallyDistributedUniformly(mesh);
}

constexpr double step_distribution_midpoint = 0.5;

double stepDistribution(const mfem::Vector& x) {
  if (x[0] > step_distribution_midpoint) {
    return 1.0;
  }
  else {
    return 0.0;
  }
}

void testThatStepDistributedPointsAllLieToTheRightOfTheMidpoint(
  mfem::Element::Type element_type,
  int num_points_to_generate = 100
) {
  auto mesh = std::make_shared<mfem::Mesh>(createMeshOfUnitBoxWith2ElemsPerDimension(element_type));

  std::default_random_engine generator;
  MeshDistribution distribution(mesh, stepDistribution);
  for (int i = 0; i < num_points_to_generate; i++) {
    mfem::Vector point;
    std::tie(point, std::ignore) = distribution.generateRandomPointAndElement(generator);
    EXPECT_GE(point[0], step_distribution_midpoint);
  }
}

TEST(MeshDistribution, TestThatStepDistributedPointsAllLieToTheRightOfTheMidpointInLineMeshes) {
  testThatStepDistributedPointsAllLieToTheRightOfTheMidpoint(mfem::Element::SEGMENT);
}

TEST(MeshDistribution, TestThatStepDistributedPointsAllLieToTheRightOfTheMidpointInQuadMeshes) {
  testThatStepDistributedPointsAllLieToTheRightOfTheMidpoint(mfem::Element::QUADRILATERAL);
}

TEST(MeshDistribution, TestThatStepDistributedPointsAllLieToTheRightOfTheMidpointInTriMeshes) {
  testThatStepDistributedPointsAllLieToTheRightOfTheMidpoint(mfem::Element::TRIANGLE);
}

TEST(MeshDistribution, TestThatStepDistributedPointsAllLieToTheRightOfTheMidpointInHexMeshes) {
  testThatStepDistributedPointsAllLieToTheRightOfTheMidpoint(mfem::Element::HEXAHEDRON);
}

TEST(MeshDistribution, TestThatStepDistributedPointsAllLieToTheRightOfTheMidpointInTetMeshes) {
  testThatStepDistributedPointsAllLieToTheRightOfTheMidpoint(mfem::Element::TETRAHEDRON);
}

constexpr double linear_distribution_slope = 2.0;

double linearDistribution(const mfem::Vector& x) {
  return linear_distribution_slope * x[0];
}

void testThatLinearlyDistributedPointsAreActuallyDistributedLinearly(
  mfem::Element::Type element_type,
  int num_elems_per_dim = 6,
  int num_points_to_generate = 20000,
  double relative_tolerance = 0.1
) {
  const int num_bins = num_elems_per_dim;
  std::vector<int> bins(num_bins, 0);

  auto mesh = std::make_shared<mfem::Mesh>();
  switch (element_type) {
  case mfem::Element::SEGMENT:
    *mesh = mfem::Mesh::MakeCartesian1D(num_elems_per_dim);
    break;
  case mfem::Element::TRIANGLE:
  case mfem::Element::QUADRILATERAL:
    *mesh = mfem::Mesh::MakeCartesian2D(num_elems_per_dim, num_elems_per_dim, element_type);
    break;
  case mfem::Element::TETRAHEDRON:
  case mfem::Element::HEXAHEDRON:
    *mesh = mfem::Mesh::MakeCartesian3D(num_elems_per_dim, num_elems_per_dim, num_elems_per_dim, element_type);
    break;
  default:
    errorWithDeveloperMessage("Element type not supported.");
  }

  std::default_random_engine generator;
  MeshDistribution distribution(mesh, linearDistribution);
  for (int i = 0; i < num_points_to_generate; i++) {
    mfem::Vector point;
    std::tie(point, std::ignore) = distribution.generateRandomPointAndElement(generator);
    bins[point[0] * num_bins] += 1;
  }

  for (int ibin = 0; ibin < num_bins; ibin++) {
    const int num_points_generated = bins[ibin];
    const double bin_left_endpoint = 1.0 * ibin / num_bins;
    const double bin_right_endpoint = (ibin + 1.0) / num_bins;
    const double expected_points_generated =
      num_points_to_generate * 0.5 * linear_distribution_slope * (std::pow(bin_right_endpoint, 2.0) - std::pow(bin_left_endpoint, 2.0));
    const double absolute_tolerance = relative_tolerance * expected_points_generated;
    EXPECT_NEAR(num_points_generated, expected_points_generated, absolute_tolerance);
  }
}

TEST(MeshDistribution, TestThatLinearlyDistributedPointsAreActuallyDistributedLinearlyInLineMesh) {
  testThatLinearlyDistributedPointsAreActuallyDistributedLinearly(mfem::Element::SEGMENT);
}

TEST(MeshDistribution, TestThatLinearlyDistributedPointsAreActuallyDistributedLinearlyInQuadMesh) {
  testThatLinearlyDistributedPointsAreActuallyDistributedLinearly(mfem::Element::QUADRILATERAL);
}

TEST(MeshDistribution, TestThatLinearlyDistributedPointsAreActuallyDistributedLinearlyInTriMesh) {
  testThatLinearlyDistributedPointsAreActuallyDistributedLinearly(mfem::Element::TRIANGLE);
}

TEST(MeshDistribution, TestThatLinearlyDistributedPointsAreActuallyDistributedLinearlyInHexMesh) {
  testThatLinearlyDistributedPointsAreActuallyDistributedLinearly(mfem::Element::HEXAHEDRON);
}

TEST(MeshDistribution, TestThatLinearlyDistributedPointsAreActuallyDistributedLinearlyInTetMesh) {
  testThatLinearlyDistributedPointsAreActuallyDistributedLinearly(mfem::Element::TETRAHEDRON);
}

} // namespace
