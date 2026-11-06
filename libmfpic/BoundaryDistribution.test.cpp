#include <libmfpic/BoundaryDistribution.hpp>
#include <libmfpic/MeshUtilities.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace {

using namespace mfpic;
using namespace testing;

struct SideOfUnitBox {
  mfem::Element::Type element_type;
  std::string side_name;
  int normal_dimension;
  double side_coordinate;
};

std::string sideErrorMessage(const SideOfUnitBox& side) {
  return "Failed with element type " + std::to_string(static_cast<int>(side.element_type)) + ", side " + side.side_name;
}

const std::vector<SideOfUnitBox> sides_to_test{
  {mfem::Element::SEGMENT, "left", 0, 0.0},
  {mfem::Element::SEGMENT, "right", 0, 1.0},
  {mfem::Element::QUADRILATERAL, "left", 0, 0.0},
  {mfem::Element::QUADRILATERAL, "top", 1, 1.0},
  {mfem::Element::TRIANGLE, "right", 0, 1.0},
  {mfem::Element::TRIANGLE, "bottom", 1, 0.0},
  {mfem::Element::HEXAHEDRON, "left", 0, 0.0},
  {mfem::Element::HEXAHEDRON, "front", 2, 1.0},
  {mfem::Element::TETRAHEDRON, "top", 1, 1.0},
  {mfem::Element::TETRAHEDRON, "back", 2, 0.0},
};

int boundaryAttributeOfSide(const mfem::Mesh& mesh, const std::string& side_name) {
  return getSideNameToBoundaryAttributeForInlineMeshes(mesh.Dimension()).at(side_name);
}

TEST(BoundaryDistribution, TotalAreaOfUnitBoxSideIsOne) {
  for (const SideOfUnitBox& side : sides_to_test) {
    SCOPED_TRACE(sideErrorMessage(side));
    auto mesh = std::make_shared<mfem::Mesh>(createMeshOfUnitBoxWith2ElemsPerDimension(side.element_type));

    const BoundaryDistribution distribution(mesh, boundaryAttributeOfSide(*mesh, side.side_name));

    EXPECT_DOUBLE_EQ(distribution.getTotalArea(), 1.0);
  }
}

TEST(BoundaryDistribution, GeneratedPointsLieOnSideWithOutwardNormalAndAdjacentElement) {
  constexpr int num_points = 200;
  constexpr double element_width = 0.5;
  for (const SideOfUnitBox& side : sides_to_test) {
    SCOPED_TRACE(sideErrorMessage(side));
    auto mesh = std::make_shared<mfem::Mesh>(createMeshOfUnitBoxWith2ElemsPerDimension(side.element_type));
    BoundaryDistribution distribution(mesh, boundaryAttributeOfSide(*mesh, side.side_name));
    RandomNumberGenerator generator;

    mfem::Vector expected_outward_unit_normal({0.0, 0.0, 0.0});
    expected_outward_unit_normal[side.normal_dimension] = side.side_coordinate == 0.0 ? -1.0 : 1.0;

    for (int i = 0; i < num_points; i++) {
      const BoundaryPoint point = distribution.generateRandomBoundaryPoint(generator);
      ASSERT_EQ(point.position.Size(), 3);
      ASSERT_EQ(point.outward_unit_normal.Size(), 3);
      EXPECT_DOUBLE_EQ(point.position[side.normal_dimension], side.side_coordinate);
      for (int dimension = 0; dimension < 3; dimension++) {
        EXPECT_DOUBLE_EQ(point.outward_unit_normal[dimension], expected_outward_unit_normal[dimension]);
      }

      mfem::Vector element_center;
      mesh->GetElementCenter(point.element, element_center);
      for (int dimension = 0; dimension < mesh->Dimension(); dimension++) {
        EXPECT_LE(std::abs(point.position[dimension] - element_center[dimension]), element_width);
      }
    }
  }
}

TEST(BoundaryDistribution, GeneratedPointsAreUniformAlongSide) {
  constexpr int num_points = 40000;
  constexpr double relative_tolerance = 0.03;
  for (const mfem::Element::Type element_type : {mfem::Element::QUADRILATERAL, mfem::Element::TRIANGLE}) {
    SCOPED_TRACE("element type " + std::to_string(static_cast<int>(element_type)));
    auto mesh = std::make_shared<mfem::Mesh>(createMeshOfUnitBoxWith2ElemsPerDimension(element_type));
    BoundaryDistribution distribution(mesh, boundaryAttributeOfSide(*mesh, "left"));
    RandomNumberGenerator generator;

    int num_points_in_lower_quarter = 0;
    double mean_y = 0.0;
    for (int i = 0; i < num_points; i++) {
      const double y = distribution.generateRandomBoundaryPoint(generator).position[1];
      mean_y += y / num_points;
      if (y < 0.25) num_points_in_lower_quarter++;
    }

    EXPECT_NEAR(mean_y, 0.5, relative_tolerance * 0.5);
    constexpr int expected_num_points_in_lower_quarter = 0.25 * num_points;
    EXPECT_NEAR(
      num_points_in_lower_quarter,
      expected_num_points_in_lower_quarter,
      relative_tolerance * expected_num_points_in_lower_quarter
    );
  }
}

TEST(BoundaryDistribution, AttributeWithNoFacesThrows) {
  auto mesh = std::make_shared<mfem::Mesh>(createMeshOfUnitBoxWith2ElemsPerDimension(mfem::Element::QUADRILATERAL));

  EXPECT_THAT(
    [&]() {BoundaryDistribution distribution(mesh, 99);},
    ThrowsMessage<std::logic_error>(HasSubstr("has no non-periodic boundary faces!"))
  );
}

} // namespace
