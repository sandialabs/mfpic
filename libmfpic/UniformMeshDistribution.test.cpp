#include <libmfpic/UniformMeshDistribution.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

void testThatRandomlyGeneratedElementsAreDistributedUniformly(
  std::shared_ptr<mfem::Mesh> mesh,
  int num_points_to_generate = 20000,
  double relative_tolerance = 0.1
) {
  const int num_elements = mesh->GetNE();
  std::vector<int> num_points_generated_in_element(num_elements, 0);

  std::mt19937 generator;
  UniformMeshDistribution distribution(mesh);
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

void testThatRandomlyGeneratedPointsAreDistributedUniformlyInASingleDimension(
  std::shared_ptr<mfem::Mesh> mesh,
  int dimension,
  int num_points_to_generate = 20000,
  int num_bins = 5,
  double relative_tolerance = 0.1
) {
  ASSERT_LT(dimension, mesh->Dimension());
  std::vector<int> bins(num_bins, 0);

  std::mt19937 generator;
  UniformMeshDistribution distribution(mesh);
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

void testThatRandomlyGeneratedPointsAreDistributedUniformly(
  std::shared_ptr<mfem::Mesh> mesh,
  int num_points_to_generate = 20000,
  int num_bins = 5,
  double relative_tolerance = 0.1
) {
  for (int dimension = 0; dimension < mesh->Dimension(); dimension++) {
    testThatRandomlyGeneratedPointsAreDistributedUniformlyInASingleDimension(
      mesh,
      dimension,
      num_points_to_generate,
      num_bins,
      relative_tolerance
    );
  }
}

TEST(UniformMeshDistribution, TestThatRandomlyGeneratedPointsAreDistributedUniformlyFor1DMeshes) {
  constexpr int num_elements = 10;
  auto mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elements));

  testThatRandomlyGeneratedElementsAreDistributedUniformly(mesh);
  testThatRandomlyGeneratedPointsAreDistributedUniformly(mesh);
}

TEST(UniformMeshDistribution, TestThatRandomlyGeneratedPointsAreDistributedUniformlyForQuadMeshes) {
  constexpr int num_elements_per_dim = 2;
  constexpr mfem::Element::Type type = mfem::Element::QUADRILATERAL;
  auto mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian2D(num_elements_per_dim, num_elements_per_dim, type));

  testThatRandomlyGeneratedElementsAreDistributedUniformly(mesh);
  testThatRandomlyGeneratedPointsAreDistributedUniformly(mesh);
}

TEST(UniformMeshDistribution, TestThatRandomlyGeneratedPointsAreDistributedUniformlyForTriMeshes) {
  constexpr int num_elements_per_dim = 2;
  constexpr mfem::Element::Type type = mfem::Element::TRIANGLE;
  auto mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian2D(num_elements_per_dim, num_elements_per_dim, type));

  testThatRandomlyGeneratedElementsAreDistributedUniformly(mesh);
  testThatRandomlyGeneratedPointsAreDistributedUniformly(mesh);
}

TEST(UniformMeshDistribution, TestThatRandomlyGeneratedPointsAreDistributedUniformlyForHexMeshes) {
  constexpr int num_elements_per_dim = 2;
  constexpr mfem::Element::Type type = mfem::Element::HEXAHEDRON;
  auto mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(num_elements_per_dim, num_elements_per_dim, num_elements_per_dim, type));

  testThatRandomlyGeneratedElementsAreDistributedUniformly(mesh);
  testThatRandomlyGeneratedPointsAreDistributedUniformly(mesh);
}

TEST(UniformMeshDistribution, TestThatRandomlyGeneratedPointsAreDistributedUniformlyForTetMeshes) {
  constexpr int num_elements_per_dim = 2;
  constexpr mfem::Element::Type type = mfem::Element::TETRAHEDRON;
  auto mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(num_elements_per_dim, num_elements_per_dim, num_elements_per_dim, type));

  testThatRandomlyGeneratedElementsAreDistributedUniformly(mesh);
  testThatRandomlyGeneratedPointsAreDistributedUniformly(mesh);
}

} // namespace
