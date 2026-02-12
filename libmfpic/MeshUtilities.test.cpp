#include <libmfpic/MeshUtilities.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

mfem::Mesh createTetMeshOfSquarePyramid() {
  constexpr int dimensions = 3;
  constexpr int initial_num_vertices = 0;
  constexpr int initial_num_elems = 0;
  constexpr int initial_num_boundary_elems = 0;
  mfem::Mesh mesh(dimensions, initial_num_vertices, initial_num_elems, initial_num_boundary_elems, dimensions);
  // mfem::Mesh mesh;
  mesh.AddVertex( 0.0,  0.0,  0.0);
  mesh.AddVertex(-1.0,  0.0,  0.0);
  mesh.AddVertex( 0.0, -1.0,  0.0);
  mesh.AddVertex( 1.0,  0.0,  0.0);
  mesh.AddVertex( 0.0,  1.0,  0.0);
  mesh.AddVertex( 0.0,  0.0,  1.0);
  mesh.AddTet(0, 1, 2, 5);
  mesh.AddTet(0, 2, 3, 5);
  mesh.AddTet(0, 3, 4, 5);
  mesh.AddTet(0, 4, 1, 5);
  mesh.GenerateBoundaryElements();
  mesh.FinalizeTetMesh();
  return mesh;
}

TEST(MeshUtilities, GetElementFaceCentroidWorksIn1D) {
  constexpr int num_elems = 4;
  mfem::Mesh non_periodic_mesh = mfem::Mesh::MakeCartesian1D(num_elems);
  std::vector<int> periodic_vertex_mapping = non_periodic_mesh.CreatePeriodicVertexMapping({mfem::Vector({1.0})});
  mfem::Mesh periodic_mesh = mfem::Mesh::MakePeriodic(non_periodic_mesh, periodic_vertex_mapping);

  constexpr double dx = 1.0 / num_elems;
  for (const mfem::Mesh* mesh : {&non_periodic_mesh, &periodic_mesh}) {
    for (int element = 0; element < num_elems; element++) {
      EXPECT_DOUBLE_EQ(getElementFaceCentroid(*mesh, element, 0)[0], element*dx);
      EXPECT_DOUBLE_EQ(getElementFaceCentroid(*mesh, element, 1)[0], (element + 1)*dx);
    }
  }
}

TEST(MeshUtilities, GetElementFaceCentroidWorksInQuadMeshes) {
  constexpr int num_x_elems = 4;
  constexpr int num_y_elems = 1;
  constexpr mfem::Element::Type element_type = mfem::Element::QUADRILATERAL;
  constexpr bool generate_edges = true;
  constexpr double domain_side_length = 1.0;
  constexpr bool space_filling_curve_ordering = false;
  mfem::Mesh non_periodic_mesh = mfem::Mesh::MakeCartesian2D(
    num_x_elems,
    num_y_elems,
    element_type,
    generate_edges,
    domain_side_length,
    domain_side_length,
    space_filling_curve_ordering
  );
  std::vector<int> periodic_vertex_mapping = non_periodic_mesh.CreatePeriodicVertexMapping({mfem::Vector({1.0, 0.0})});
  mfem::Mesh periodic_mesh = mfem::Mesh::MakePeriodic(non_periodic_mesh, periodic_vertex_mapping);

  // Quad face ordering
  //    2
  //  *---*
  // 3|   |1
  //  *---*
  //    0
  constexpr double dx = domain_side_length / num_x_elems;
  for (const mfem::Mesh* mesh : {&non_periodic_mesh, &periodic_mesh}) {
    for (int element = 0; element < num_x_elems; element++) {
      const double element_center_x = (element + 0.5) * dx;
      EXPECT_DOUBLE_EQ(getElementFaceCentroid(*mesh, element, 0)[0], element_center_x);
      EXPECT_DOUBLE_EQ(getElementFaceCentroid(*mesh, element, 0)[1], 0.0);
      EXPECT_DOUBLE_EQ(getElementFaceCentroid(*mesh, element, 1)[0], (element + 1)*dx);
      EXPECT_DOUBLE_EQ(getElementFaceCentroid(*mesh, element, 1)[1], domain_side_length / 2.0);
      EXPECT_DOUBLE_EQ(getElementFaceCentroid(*mesh, element, 2)[0], element_center_x);
      EXPECT_DOUBLE_EQ(getElementFaceCentroid(*mesh, element, 2)[1], domain_side_length);
      EXPECT_DOUBLE_EQ(getElementFaceCentroid(*mesh, element, 3)[0], element * dx);
      EXPECT_DOUBLE_EQ(getElementFaceCentroid(*mesh, element, 3)[1], domain_side_length / 2.0);
    }
  }
}

TEST(MeshUtilities, GetElementFaceCentroidWorksInTriMeshes) {
  constexpr int num_elems_per_dim = 1;
  constexpr mfem::Element::Type element_type = mfem::Element::TRIANGLE;
  constexpr bool generate_edges = true;
  constexpr double domain_side_length = 1.0;
  constexpr bool space_filling_curve_ordering = false;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
    num_elems_per_dim,
    num_elems_per_dim,
    element_type,
    generate_edges,
    domain_side_length,
    domain_side_length,
    space_filling_curve_ordering
  );


  // Mesh element ordering
  // *---*
  // |0/1|
  // *---*
  // Tri face ordering
  //    1
  //  *--*    0 *
  // 2| /     / |2
  //  *  0   *--*
  //           1
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 0)[0], domain_side_length / 2.0);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 0)[1], domain_side_length / 2.0);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 1)[0], domain_side_length / 2.0);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 1)[1], domain_side_length);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 2)[0], 0.0);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 2)[1], domain_side_length / 2.0);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 1, 0)[0], domain_side_length / 2.0);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 1, 0)[1], domain_side_length / 2.0);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 1, 1)[0], domain_side_length / 2.0);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 1, 1)[1], 0.0);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 1, 2)[0], domain_side_length);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 1, 2)[1], domain_side_length / 2.0);
}

TEST(MeshUtilities, GetElementFaceCentroidWorksInHexMeshes) {
  constexpr int num_elems_per_dim = 1;
  constexpr mfem::Element::Type element_type = mfem::Element::HEXAHEDRON;
  constexpr double domain_side_length = 1.0;
  constexpr bool space_filling_curve_ordering = false;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
    num_elems_per_dim,
    num_elems_per_dim,
    num_elems_per_dim,
    element_type,
    domain_side_length,
    domain_side_length,
    domain_side_length,
    space_filling_curve_ordering
  );

  // Element face ordering
  // z
  // |_ x
  //     5
  //   *---*
  // 4 |   | 2
  //   *---*
  //     0
  // y
  // |_ x
  //     3
  //   *---*
  // 4 |   | 2
  //   *---*
  //     1
  constexpr double side_center = domain_side_length / 2.0;
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 0)[0], side_center);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 0)[1], side_center);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 0)[2], 0.0);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 1)[0], side_center);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 1)[1], 0.0);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 1)[2], side_center);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 2)[0], domain_side_length);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 2)[1], side_center);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 2)[2], side_center);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 3)[0], side_center);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 3)[1], domain_side_length);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 3)[2], side_center);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 4)[0], 0.0);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 4)[1], side_center);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 4)[2], side_center);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 5)[0], side_center);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 5)[1], side_center);
  EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, 0, 5)[2], domain_side_length);
}

TEST(MeshUtilities, GetElementFaceCentroidWorksInTetMeshes) {
  mfem::Mesh mesh = createTetMeshOfSquarePyramid();

  constexpr double one_third = 1.0 / 3.0;
  for (int element = 0; element < mesh.GetNE(); element++) {
    mfem::Vector element_center;
    mesh.GetElementCenter(element, element_center);
    mfem::Vector bottom_face_centroid{
      std::copysign(one_third, element_center[0]),
      std::copysign(one_third, element_center[1]),
      0.0
    };
    mfem::Vector top_face_centroid{
      std::copysign(one_third, element_center[0]),
      std::copysign(one_third, element_center[1]),
      one_third
    };
    mfem::Vector vertical_face_centroid{
      0.0,
      std::copysign(one_third, element_center[1]),
      one_third
    };
    mfem::Vector horizontal_face_centroid{
      std::copysign(one_third, element_center[0]),
      0.0,
      one_third
    };
    const bool face_1_is_vertical = element == 0 or element == 2;
    mfem::Vector& face_1_centroid = face_1_is_vertical ? vertical_face_centroid : horizontal_face_centroid;
    mfem::Vector& face_2_centroid = face_1_is_vertical ? horizontal_face_centroid : vertical_face_centroid;
    for (int idim = 0; idim < 3; idim++) {
      EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, element, 0)[idim], top_face_centroid[idim]);
      EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, element, 1)[idim], face_1_centroid[idim]);
      EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, element, 2)[idim], face_2_centroid[idim]);
      EXPECT_DOUBLE_EQ(getElementFaceCentroid(mesh, element, 3)[idim], bottom_face_centroid[idim]);
    }
  }
}

TEST(MeshUtilities, GetElementOnOtherSideOfFaceWorksIn1D) {
  constexpr int num_elems = 4;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);

  for (int element = 0; element < num_elems - 1; element++) {
    EXPECT_EQ(getElementOnOtherSideOfFace(mesh, element, 1), element+1);
    EXPECT_EQ(getElementOnOtherSideOfFace(mesh, element+1, 0), element);
  }
}

TEST(MeshUtilities, ElementOnOtherSideOfBoundaryElementsAreInvalidIn1D) {
  constexpr int num_elems = 4;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);

  EXPECT_LT(getElementOnOtherSideOfFace(mesh, 0, 0), 0);
  EXPECT_LT(getElementOnOtherSideOfFace(mesh, num_elems - 1, 1), 0);
}

TEST(MeshUtilities, GetElementOnOtherSideOfFaceWorksForQuadMesh) {
  constexpr int num_elems_per_dim = 3;
  constexpr bool generate_edges = true;
  constexpr mfem::Element::Type element_type = mfem::Element::QUADRILATERAL;
  constexpr double domain_side_length = 1.0;
  constexpr bool space_filling_curve_ordering = false;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
    num_elems_per_dim,
    num_elems_per_dim,
    element_type,
    generate_edges,
    domain_side_length,
    domain_side_length,
    space_filling_curve_ordering
  );

  // Mesh element ordering
  // *---*---*---*
  // | 6 | 7 | 8 |
  // *---*---*---*
  // | 3 | 4 | 5 |
  // *---*---*---*
  // | 0 | 1 | 2 |
  // *---*---*---*

  // Quad face ordering
  //    2
  //  *---*
  // 3|   |1
  //  *---*
  //    0

  // Corners
  EXPECT_LT(getElementOnOtherSideOfFace(mesh, 0, 0), 0);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 0, 1), 1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 0, 2), 3);
  EXPECT_LT(getElementOnOtherSideOfFace(mesh, 0, 3), 0);
  EXPECT_LT(getElementOnOtherSideOfFace(mesh, 2, 0), 0);
  EXPECT_LT(getElementOnOtherSideOfFace(mesh, 2, 1), 0);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 2, 2), 5);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 2, 3), 1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 6, 0), 3);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 6, 1), 7);
  EXPECT_LT(getElementOnOtherSideOfFace(mesh, 6, 2), 0);
  EXPECT_LT(getElementOnOtherSideOfFace(mesh, 6, 3), 0);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 8, 0), 5);
  EXPECT_LT(getElementOnOtherSideOfFace(mesh, 8, 1), 0);
  EXPECT_LT(getElementOnOtherSideOfFace(mesh, 8, 2), 0);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 8, 3), 7);
  // Edges
  EXPECT_LT(getElementOnOtherSideOfFace(mesh, 1, 0), 0);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 1, 1), 2);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 1, 2), 4);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 1, 3), 0);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 3, 0), 0);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 3, 1), 4);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 3, 2), 6);
  EXPECT_LT(getElementOnOtherSideOfFace(mesh, 3, 3), 0);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 5, 0), 2);
  EXPECT_LT(getElementOnOtherSideOfFace(mesh, 5, 1), 0);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 5, 2), 8);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 5, 3), 4);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 7, 0), 4);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 7, 1), 8);
  EXPECT_LT(getElementOnOtherSideOfFace(mesh, 7, 2), 0);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 7, 3), 6);
  // Center
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 4, 0), 1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 4, 1), 5);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 4, 2), 7);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 4, 3), 3);
}

TEST(MeshUtilities, GetElementOnOtherSideOfFaceWorksForTriMesh) {
  constexpr int num_elems_per_dim = 2;
  constexpr bool generate_edges = true;
  constexpr mfem::Element::Type element_type = mfem::Element::TRIANGLE;
  constexpr double domain_side_length = 1.0;
  constexpr bool space_filling_curve_ordering = false;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
    num_elems_per_dim,
    num_elems_per_dim,
    element_type,
    generate_edges,
    domain_side_length,
    domain_side_length,
    space_filling_curve_ordering
  );

  // Mesh element ordering
  // *---*---*
  // |4/5|6/7|
  // *---*---*
  // |0/1|2/3|
  // *---*---*

  // Tri face ordering
  //    1
  //  *--*    0 *
  // 2| /     / |2
  //  *  0   *--*
  //           1

  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 0, 0), 1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 0, 1), 5);
  EXPECT_LT(getElementOnOtherSideOfFace(mesh, 0, 2), 0);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 1, 0), 0);
  EXPECT_LT(getElementOnOtherSideOfFace(mesh, 1, 1), 0);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 1, 2), 2);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 2, 0), 3);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 2, 1), 7);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 2, 2), 1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 3, 0), 2);
  EXPECT_LT(getElementOnOtherSideOfFace(mesh, 3, 1), 0);
  EXPECT_LT(getElementOnOtherSideOfFace(mesh, 3, 2), 0);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 4, 0), 5);
  EXPECT_LT(getElementOnOtherSideOfFace(mesh, 4, 1), 0);
  EXPECT_LT(getElementOnOtherSideOfFace(mesh, 4, 2), 0);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 5, 0), 4);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 5, 1), 0);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 5, 2), 6);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 6, 0), 7);
  EXPECT_LT(getElementOnOtherSideOfFace(mesh, 6, 1), 0);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 6, 2), 5);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 7, 0), 6);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 7, 1), 2);
  EXPECT_LT(getElementOnOtherSideOfFace(mesh, 7, 2), 0);
}

TEST(MeshUtilities, GetElementOnOtherSideOfFaceWorksForHexMesh) {
  constexpr int num_elems_per_dim = 2;
  constexpr mfem::Element::Type element_type = mfem::Element::HEXAHEDRON;
  constexpr double domain_side_length = 1.0;
  constexpr bool space_filling_curve_ordering = false;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
    num_elems_per_dim,
    num_elems_per_dim,
    num_elems_per_dim,
    element_type,
    domain_side_length,
    domain_side_length,
    domain_side_length,
    space_filling_curve_ordering
  );

  // Mesh element ordering
  // Layer 0
  // *---*---*
  // | 2 | 3 |
  // *---*---*
  // | 0 | 1 |
  // *---*---*
  // Layer 1
  // *---*---*
  // | 6 | 7 |
  // *---*---*
  // | 4 | 5 |
  // *---*---*

  // Element face ordering
  // z
  // |_ x
  //     5
  //   *---*
  // 4 |   | 2
  //   *---*
  //     0
  // y
  // |_ x
  //     3
  //   *---*
  // 4 |   | 2
  //   *---*
  //     1

  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 0, 0), -1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 0, 1), -1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 0, 2),  1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 0, 3),  2);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 0, 4), -1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 0, 5),  4);

  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 1, 0), -1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 1, 1), -1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 1, 2), -1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 1, 3),  3);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 1, 4),  0);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 1, 5),  5);

  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 2, 0), -1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 2, 1),  0);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 2, 2),  3);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 2, 3), -1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 2, 4), -1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 2, 5),  6);

  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 3, 0), -1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 3, 1),  1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 3, 2), -1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 3, 3), -1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 3, 4),  2);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 3, 5),  7);

  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 4, 0),  0);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 4, 1), -1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 4, 2),  5);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 4, 3),  6);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 4, 4), -1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 4, 5), -1);

  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 5, 0),  1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 5, 1), -1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 5, 2), -1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 5, 3),  7);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 5, 4),  4);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 5, 5), -1);

  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 6, 0),  2);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 6, 1),  4);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 6, 2),  7);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 6, 3), -1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 6, 4), -1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 6, 5), -1);

  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 7, 0),  3);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 7, 1),  5);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 7, 2), -1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 7, 3), -1);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 7, 4),  6);
  EXPECT_EQ(getElementOnOtherSideOfFace(mesh, 7, 5), -1);
}

TEST(MeshUtilities, GetElementOnOtherSideOfFaceWorksForTetMesh) {
  mfem::Mesh mesh = createTetMeshOfSquarePyramid();

  const int num_elems = mesh.GetNE();
  EXPECT_EQ(num_elems, 4);
  for (int element = 0; element < mesh.GetNE(); element++) {
    EXPECT_EQ(getElementOnOtherSideOfFace(mesh, element, 0), -1);
    EXPECT_EQ(getElementOnOtherSideOfFace(mesh, element, 1), (element + 1) % num_elems);
    EXPECT_EQ(getElementOnOtherSideOfFace(mesh, element, 2), (element + num_elems - 1) % num_elems);
    EXPECT_EQ(getElementOnOtherSideOfFace(mesh, element, 3), -1);
  }
}

TEST(MeshUtilities, OutwardFacingUnitNormalsComputedCorrectlyIn1D) {
  constexpr int num_elems = 4;
  mfem::Mesh non_periodic_mesh = mfem::Mesh::MakeCartesian1D(num_elems);
  std::vector<int> periodic_vertex_mapping = non_periodic_mesh.CreatePeriodicVertexMapping({mfem::Vector({1.0})});
  mfem::Mesh periodic_mesh = mfem::Mesh::MakePeriodic(non_periodic_mesh, periodic_vertex_mapping);

  for (const mfem::Mesh* mesh : {&non_periodic_mesh, &periodic_mesh}) {
    for (int element = 0; element < num_elems; element++) {
      EXPECT_DOUBLE_EQ(getElementFaceOutwardUnitNormal(*mesh, element, 0)[0], -1.0);
      EXPECT_DOUBLE_EQ(getElementFaceOutwardUnitNormal(*mesh, element, 1)[0],  1.0);
    }
  }
}

TEST(MeshUtilities, OutwardFacingUnitNormalsComputedCorrectlyForQuads) {
  constexpr int num_elems_per_dim = 3;
  constexpr bool generate_edges = true;
  constexpr mfem::Element::Type element_type = mfem::Element::QUADRILATERAL;
  mfem::Mesh non_periodic_mesh = mfem::Mesh::MakeCartesian2D(
    num_elems_per_dim,
    num_elems_per_dim,
    element_type,
    generate_edges
  );
  std::vector<int> periodic_vertex_mapping = non_periodic_mesh.CreatePeriodicVertexMapping({
    mfem::Vector({1.0, 0.0}),
    mfem::Vector({0.0, 1.0})
  });
  mfem::Mesh periodic_mesh = mfem::Mesh::MakePeriodic(non_periodic_mesh, periodic_vertex_mapping);

  const std::vector<mfem::Vector> expected_unit_normals = {
    mfem::Vector({ 0.0, -1.0}),
    mfem::Vector({ 1.0,  0.0}),
    mfem::Vector({ 0.0,  1.0}),
    mfem::Vector({-1.0,  0.0})
  };
  for (const mfem::Mesh* mesh : {&non_periodic_mesh, &periodic_mesh}) {
    for (int element = 0; element < mesh->GetNE(); element++) {
      for (int face = 0; face < 4; face++) {
        for (int dimension = 0; dimension < 2; dimension++) {
          EXPECT_DOUBLE_EQ(
            getElementFaceOutwardUnitNormal(*mesh, element, face)[dimension],
            expected_unit_normals[face][dimension]
          );
        }
      }
    }
  }
}

TEST(MeshUtilities, OutwardFacingUnitNormalsComputedCorrectlyForTris) {
  constexpr int num_elems_per_dim = 3;
  constexpr bool generate_edges = true;
  constexpr mfem::Element::Type element_type = mfem::Element::TRIANGLE;
  mfem::Mesh non_periodic_mesh = mfem::Mesh::MakeCartesian2D(
    num_elems_per_dim,
    num_elems_per_dim,
    element_type,
    generate_edges
  );
  std::vector<int> periodic_vertex_mapping = non_periodic_mesh.CreatePeriodicVertexMapping({
    mfem::Vector({1.0, 0.0}),
    mfem::Vector({0.0, 1.0})
  });
  mfem::Mesh periodic_mesh = mfem::Mesh::MakePeriodic(non_periodic_mesh, periodic_vertex_mapping);

  const std::vector<mfem::Vector> even_element_expected_unit_normals = {
    mfem::Vector({M_SQRT1_2, -M_SQRT1_2}),
    mfem::Vector({0.0, 1.0}),
    mfem::Vector({-1.0, 0.0})
  };
  const std::vector<mfem::Vector> odd_element_expected_unit_normals = {
    mfem::Vector({-M_SQRT1_2, M_SQRT1_2}),
    mfem::Vector({0.0, -1.0}),
    mfem::Vector({1.0, 0.0})
  };
  for (const mfem::Mesh* mesh : {&non_periodic_mesh, &periodic_mesh}) {
    for (int element = 0; element < mesh->GetNE(); element++) {
      for (int face = 0; face < 3; face++) {
        for (int dimension = 0; dimension < 2; dimension++) {
          if (element % 2) {
            EXPECT_DOUBLE_EQ(
              getElementFaceOutwardUnitNormal(*mesh, element, face)[dimension],
              odd_element_expected_unit_normals[face][dimension]
            );
          } else {
            EXPECT_DOUBLE_EQ(
              getElementFaceOutwardUnitNormal(*mesh, element, face)[dimension],
              even_element_expected_unit_normals[face][dimension]
            );
          }
        }
      }
    }
  }
}

TEST(MeshUtilities, OutwardFacingUnitNormalsComputedCorrectlyForHexes) {
  constexpr int num_elems_per_dim = 3;
  constexpr mfem::Element::Type element_type = mfem::Element::HEXAHEDRON;
  mfem::Mesh non_periodic_mesh = mfem::Mesh::MakeCartesian3D(
    num_elems_per_dim,
    num_elems_per_dim,
    num_elems_per_dim,
    element_type
  );
  std::vector<int> periodic_vertex_mapping = non_periodic_mesh.CreatePeriodicVertexMapping({
    mfem::Vector({1.0, 0.0, 0.0}),
    mfem::Vector({0.0, 1.0, 0.0}),
    mfem::Vector({0.0, 0.0, 1.0})
  });
  mfem::Mesh periodic_mesh = mfem::Mesh::MakePeriodic(non_periodic_mesh, periodic_vertex_mapping);

  const std::vector<mfem::Vector> expected_unit_normals = {
    mfem::Vector({ 0.0,  0.0, -1.0}),
    mfem::Vector({ 0.0, -1.0,  0.0}),
    mfem::Vector({ 1.0,  0.0,  0.0}),
    mfem::Vector({ 0.0,  1.0,  0.0}),
    mfem::Vector({-1.0,  0.0,  0.0}),
    mfem::Vector({ 0.0,  0.0,  1.0})
  };
  for (const mfem::Mesh* mesh : {&non_periodic_mesh, &periodic_mesh}) {
    for (int element = 0; element < mesh->GetNE(); element++) {
      for (int face = 0; face < 6; face++) {
        for (int dimension = 0; dimension < 3; dimension++) {
          EXPECT_DOUBLE_EQ(
            getElementFaceOutwardUnitNormal(*mesh, element, face)[dimension],
            expected_unit_normals[face][dimension]
          );
        }
      }
    }
  }
}

TEST(MeshUtilities, OutwardFacingUnitNormalsComputedCorrectlyForTets) {
  mfem::Mesh mesh = createTetMeshOfSquarePyramid();

  const int num_elems = mesh.GetNE();
  EXPECT_EQ(num_elems, 4);
  mfem::Vector downward_unit_normal{0.0, 0.0, -1.0};
  for (int element = 0; element < mesh.GetNE(); element++) {
    mfem::Vector element_center;
    mesh.GetElementCenter(element, element_center);
    mfem::Vector next_element_center;
    mesh.GetElementCenter((element + 1) % num_elems, next_element_center);
    mfem::Vector prev_element_center;
    mesh.GetElementCenter((element + num_elems - 1) % num_elems, prev_element_center);
    mfem::Vector face_1_outward_unit_normal = next_element_center;
    face_1_outward_unit_normal.Add(-1.0, element_center);
    face_1_outward_unit_normal /= face_1_outward_unit_normal.Norml2();
    mfem::Vector face_2_outward_unit_normal = prev_element_center;
    face_2_outward_unit_normal.Add(-1.0, element_center);
    face_2_outward_unit_normal /= face_2_outward_unit_normal.Norml2();
    mfem::Vector upward_face_unit_normal = element_center;
    upward_face_unit_normal[2] = std::abs(upward_face_unit_normal[0]);
    upward_face_unit_normal /= upward_face_unit_normal.Norml2();

    for (int idim = 0; idim < 3; idim++) {
      EXPECT_DOUBLE_EQ(getElementFaceOutwardUnitNormal(mesh, element, 0)[idim], upward_face_unit_normal[idim]);
      EXPECT_DOUBLE_EQ(getElementFaceOutwardUnitNormal(mesh, element, 1)[idim], face_1_outward_unit_normal[idim]);
      EXPECT_DOUBLE_EQ(getElementFaceOutwardUnitNormal(mesh, element, 2)[idim], face_2_outward_unit_normal[idim]);
      EXPECT_DOUBLE_EQ(getElementFaceOutwardUnitNormal(mesh, element, 3)[idim], downward_unit_normal[idim]);
    }
  }
}

void testThatGetElementFaceOfBoundaryElementReturnsUniquePairThatIsActuallyOnBoundary(const mfem::Mesh& mesh) {
  std::set<std::tuple<int, int, bool>> results_already_returned;

  for (int boundary_element = 0; boundary_element < mesh.GetNBE(); boundary_element++) {
    std::tuple<int, int, bool> element_face_result = getElementFaceOfBoundaryElement(mesh, boundary_element);
    const auto insert_result = results_already_returned.insert(element_face_result);
    const bool result_did_not_previously_exist = insert_result.second;
    EXPECT_TRUE(result_did_not_previously_exist);
    EXPECT_TRUE(std::get<2>(element_face_result));
    EXPECT_EQ(getElementOnOtherSideOfFace(mesh, std::get<0>(element_face_result), std::get<1>(element_face_result)), -1);
  }
}

TEST(MeshUtilities, Each1DMeshBoundaryElementMapsToUniqueElementFacePairThatIsOnBoundary) {
  constexpr int num_elems = 4;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);

  ASSERT_EQ(mesh.GetNBE(), 2);
  testThatGetElementFaceOfBoundaryElementReturnsUniquePairThatIsActuallyOnBoundary(mesh);
}

TEST(MeshUtilities, EachQuadMeshBoundaryElementMapsToUniqueElementFacePairThatIsOnBoundary) {
  constexpr int num_elems_per_dim = 3;
  constexpr bool generate_edges = true;
  constexpr mfem::Element::Type element_type = mfem::Element::QUADRILATERAL;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
    num_elems_per_dim,
    num_elems_per_dim,
    element_type,
    generate_edges
  );

  ASSERT_EQ(mesh.GetNBE(), 4 * num_elems_per_dim);
  testThatGetElementFaceOfBoundaryElementReturnsUniquePairThatIsActuallyOnBoundary(mesh);
}

TEST(MeshUtilities, EachTriMeshBoundaryElementMapsToUniqueElementFacePairThatIsOnBoundary) {
  constexpr int num_elems_per_dim = 3;
  constexpr bool generate_edges = true;
  constexpr mfem::Element::Type element_type = mfem::Element::TRIANGLE;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
    num_elems_per_dim,
    num_elems_per_dim,
    element_type,
    generate_edges
  );

  ASSERT_EQ(mesh.GetNBE(), 4 * num_elems_per_dim);
  testThatGetElementFaceOfBoundaryElementReturnsUniquePairThatIsActuallyOnBoundary(mesh);
}

TEST(MeshUtilities, EachHexMeshBoundaryElementMapsToUniqueElementFacePairThatIsOnBoundary) {
  constexpr int num_elems_per_dim = 3;
  constexpr mfem::Element::Type element_type = mfem::Element::HEXAHEDRON;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
    num_elems_per_dim,
    num_elems_per_dim,
    num_elems_per_dim,
    element_type
  );

  ASSERT_EQ(mesh.GetNBE(), 6 * num_elems_per_dim * num_elems_per_dim);
  testThatGetElementFaceOfBoundaryElementReturnsUniquePairThatIsActuallyOnBoundary(mesh);
}

TEST(MeshUtilities, EachTetMeshBoundaryElementMapsToUniqueElementFacePairThatIsOnBoundary) {
  mfem::Mesh mesh = createTetMeshOfSquarePyramid();

  ASSERT_EQ(mesh.GetNBE(), 8);
  testThatGetElementFaceOfBoundaryElementReturnsUniquePairThatIsActuallyOnBoundary(mesh);
}

TEST(MeshUtilities, GetMeshVolumeWorksFor1DMeshes) {
  constexpr int num_elems = 4;
  constexpr double domain_length = 123.0;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems, domain_length);

  ASSERT_DOUBLE_EQ(domain_length, getMeshVolume(mesh));
}

TEST(MeshUtilities, GetMeshVolumeWorksFor2DMeshes) {
  constexpr int num_elems_per_dim = 4;
  constexpr bool generate_edges = true;
  constexpr double x_side_length = 123.0;
  constexpr double y_side_length = 235.0;
  for (mfem::Element::Type element_type : {mfem::Element::QUADRILATERAL, mfem::Element::TRIANGLE}) {
    mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
      num_elems_per_dim,
      num_elems_per_dim,
      element_type,
      generate_edges,
      x_side_length,
      y_side_length
    );

    EXPECT_DOUBLE_EQ(x_side_length * y_side_length, getMeshVolume(mesh));
  }
}

TEST(MeshUtilities, GetMeshVolumeWorksFor3DMeshes) {
  constexpr int num_elems_per_dim = 4;
  constexpr double x_side_length = 123.0;
  constexpr double y_side_length = 235.0;
  constexpr double z_side_length = 14.4;
  for (mfem::Element::Type element_type : {mfem::Element::HEXAHEDRON, mfem::Element::TETRAHEDRON}) {
    mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
      num_elems_per_dim,
      num_elems_per_dim,
      num_elems_per_dim,
      element_type,
      x_side_length,
      y_side_length,
      z_side_length
    );

    EXPECT_DOUBLE_EQ(x_side_length * y_side_length * z_side_length, getMeshVolume(mesh));
  }
}

TEST(MeshUtilities, GetNumFacesOnElementWorksFor1DMeshes) {
  constexpr int num_elems = 4;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);

  for (int element = 0; element < mesh.GetNE(); element++) {
    EXPECT_EQ(getNumFacesOnElement(mesh, element), 2);
  }
}

TEST(MeshUtilities, GetNumFacesOnElementWorksFor2DMeshes) {
  constexpr int num_elems_per_dim = 4;
  constexpr bool generate_edges = true;
  constexpr double x_side_length = 123.0;
  constexpr double y_side_length = 235.0;
  const mfem::Element::Type element_types[] = {mfem::Element::QUADRILATERAL, mfem::Element::TRIANGLE};
  const int expected_num_faces[] = {4, 3};
  for (int i = 0; i < 2; i++) {
    mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
      num_elems_per_dim,
      num_elems_per_dim,
      element_types[i],
      generate_edges,
      x_side_length,
      y_side_length
    );

    for (int element = 0; element < mesh.GetNE(); element++) {
      EXPECT_EQ(getNumFacesOnElement(mesh, element), expected_num_faces[i]);
    }
  }
}

TEST(MeshUtilities, GetNumFacesOnElementWorksFor3DMeshes) {
  constexpr int num_elems_per_dim = 4;
  constexpr double x_side_length = 123.0;
  constexpr double y_side_length = 235.0;
  constexpr double z_side_length = 14.4;
  const mfem::Element::Type element_types[] = {mfem::Element::HEXAHEDRON, mfem::Element::TETRAHEDRON};
  const int expected_num_faces[] = {6, 4};
  for (int i = 0; i < 2; i++) {
    mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
      num_elems_per_dim,
      num_elems_per_dim,
      num_elems_per_dim,
      element_types[i],
      x_side_length,
      y_side_length,
      z_side_length
    );

    for (int element = 0; element < mesh.GetNE(); element++) {
      EXPECT_EQ(getNumFacesOnElement(mesh, element), expected_num_faces[i]);
    }
  }
}

void checkSideNameToBoundaryAttributeMap(
  const mfem::Mesh& mesh,
  const std::vector<std::string>& expected_side_names,
  const double length)
{
  std::unordered_map<std::string, int> side_name_to_boundary_attribute = getSideNameToBoundaryAttributeForInlineMeshes(
    mesh.Dimension());

  EXPECT_EQ(side_name_to_boundary_attribute.size(), expected_side_names.size());
  for (const std::string& expected_side_name : expected_side_names) {
    EXPECT_TRUE(side_name_to_boundary_attribute.contains(expected_side_name));
  }

  const int num_boundary_elements = mesh.GetNBE();
  for (const auto& [side_name, boundary_attribute] : side_name_to_boundary_attribute) {
    for (int i_boundary_element = 0; i_boundary_element < num_boundary_elements; ++i_boundary_element) {
      const mfem::Element* boundary_element = mesh.GetBdrElement(i_boundary_element);
      const int i_boundary_attribute = boundary_element->GetAttribute();
      if (i_boundary_attribute == boundary_attribute) {
        mfem::Array<int> vertex_indices;
        boundary_element->GetVertices(vertex_indices);
        for (const int& i_vertex : vertex_indices) {
          const double* vertex_position = mesh.GetVertex(i_vertex);
          if (side_name == "left") {
            EXPECT_EQ(vertex_position[0], 0.);
          } else if (side_name == "right") {
            EXPECT_EQ(vertex_position[0], length);
          } else if (side_name == "bottom") {
            EXPECT_EQ(vertex_position[1], 0.);
          } else if (side_name == "top") {
            EXPECT_EQ(vertex_position[1], length);
          } else if (side_name == "back") {
            EXPECT_EQ(vertex_position[2], 0.);
          } else if (side_name == "front") {
            EXPECT_EQ(vertex_position[2], length);
          }
        }
      }
    }
  }
}

TEST(MeshUtilities, SideStringsToAttributes1D) {
  constexpr int num_elems = 10;
  constexpr double length = 1.2;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems, length);

  checkSideNameToBoundaryAttributeMap(mesh, {"left", "right"}, length);
}

TEST(MeshUtilities, SideStringsToAttributes2DQuad) {
  constexpr int num_elems_per_dim = 3;
  constexpr mfem::Element::Type element_type = mfem::Element::QUADRILATERAL;
  constexpr bool generate_edges = true;
  constexpr double length = 2.8;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
    num_elems_per_dim,
    num_elems_per_dim,
    element_type,
    generate_edges,
    length,
    length);

  checkSideNameToBoundaryAttributeMap(mesh, {"left", "right", "bottom", "top"}, length);
}

TEST(MeshUtilities, SideStringsToAttributes2DTri) {
  constexpr int num_elems_per_dim = 4;
  constexpr mfem::Element::Type element_type = mfem::Element::TRIANGLE;
  constexpr bool generate_edges = true;
  constexpr double length = 0.8;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
    num_elems_per_dim,
    num_elems_per_dim,
    element_type,
    generate_edges,
    length,
    length);

  checkSideNameToBoundaryAttributeMap(mesh, {"left", "right", "bottom", "top"}, length);
}

TEST(MeshUtilities, SideStringsToAttributes3DHexahedron) {
  constexpr int num_elems_per_dim = 3;
  constexpr mfem::Element::Type element_type = mfem::Element::HEXAHEDRON;
  constexpr double length = 0.7;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
    num_elems_per_dim,
    num_elems_per_dim,
    num_elems_per_dim,
    element_type,
    length,
    length,
    length);

  checkSideNameToBoundaryAttributeMap(mesh, {"left", "right", "bottom", "top", "front", "back"}, length);
}

TEST(MeshUtilities, SideStringsToAttributes3DTetrahedron) {
  constexpr int num_elems_per_dim = 5;
  constexpr mfem::Element::Type element_type = mfem::Element::TETRAHEDRON;
  constexpr double length = 0.6;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
    num_elems_per_dim,
    num_elems_per_dim,
    num_elems_per_dim,
    element_type,
    length,
    length,
    length);

  checkSideNameToBoundaryAttributeMap(mesh, {"left", "right", "bottom", "top", "front", "back"}, length);
}

} // namespace
