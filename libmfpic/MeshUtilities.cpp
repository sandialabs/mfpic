#include <libmfpic/Errors.hpp>
#include <libmfpic/MeshUtilities.hpp>
#include <limits>

namespace mfpic {

namespace {

mfem::Array<int> elementFaceToMeshFaceMap(const mfem::Mesh& mesh, int element) {
  assert(element < mesh.GetNE());

  mfem::Array<int> element_face_to_mesh_face;
  if (mesh.Dimension() == 1) {
    mesh.GetElementVertices(element, element_face_to_mesh_face);
  } else if (mesh.Dimension() == 2) {
    mfem::Array<int> orientations;
    mesh.GetElementEdges(element, element_face_to_mesh_face, orientations);
  } else if (mesh.Dimension() == 3) {
    mfem::Array<int> orientations;
    mesh.GetElementFaces(element, element_face_to_mesh_face, orientations);
  } else {
    std::string error_message = "Mesh dimension must be between 1 and 3 (inclusive)!";
    errorWithDeveloperMessage(error_message);
  }

  return element_face_to_mesh_face;
}

} // namespace

mfem::Vector getElementFaceCentroid(const mfem::Mesh& mesh, int element, int element_face) {
  mfem::Array<int> element_face_to_mesh_face = elementFaceToMeshFaceMap(mesh, element);
  assert(element_face < element_face_to_mesh_face.Size());

  const int mesh_face = element_face_to_mesh_face[element_face];
  mfem::FaceElementTransformations face_element_transformations;
  mfem::IsoparametricTransformation inner_element_transformation, outer_element_transformation;
  mesh.GetFaceElementTransformations(
    mesh_face,
    face_element_transformations,
    inner_element_transformation,
    outer_element_transformation
  );
  const bool this_element_is_inner_element = element == face_element_transformations.Elem1No;
  mfem::IsoparametricTransformation& element_transformation = this_element_is_inner_element ? inner_element_transformation : outer_element_transformation;
  mfem::IntegrationPointTransformation& face_to_element_transformation = this_element_is_inner_element ? face_element_transformations.GetIntPoint1Transformation() : face_element_transformations.GetIntPoint2Transformation();

  mfem::IntegrationPoint face_centroid_in_element_reference_coords;
  mfem::Geometry::Type face_geometry = mesh.GetFaceGeometry(mesh_face);
  face_to_element_transformation.Transform(mfem::Geometries.GetCenter(face_geometry), face_centroid_in_element_reference_coords);

  mfem::Vector face_centroid;
  element_transformation.Transform(face_centroid_in_element_reference_coords, face_centroid);

  assert(face_centroid.Size() == mesh.Dimension());

  return face_centroid;
}

int getElementOnOtherSideOfFace(const mfem::Mesh& mesh, int element, int element_face) {
  mfem::Array<int> element_face_to_mesh_face = elementFaceToMeshFaceMap(mesh, element);
  assert(element_face < element_face_to_mesh_face.Size());

  int inner_element, outer_element;
  mesh.GetFaceElements(element_face_to_mesh_face[element_face], &inner_element, &outer_element);

  return element == inner_element ? outer_element : inner_element;
}

mfem::Vector getElementFaceOutwardUnitNormal(const mfem::Mesh& mesh, int element, int element_face) {
  mfem::Array<int> element_face_to_mesh_face = elementFaceToMeshFaceMap(mesh, element);
  assert(element_face < element_face_to_mesh_face.Size());
  mfem::Vector unit_normal(mesh.Dimension());
  if (mesh.Dimension() == 1) {
    mfem::Vector element_face_centroid = getElementFaceCentroid(mesh, element, element_face);
    mfem::Vector element_other_face_centroid = getElementFaceCentroid(mesh, element, 1 - element_face);
    unit_normal[0] = element_face_centroid[0] > element_other_face_centroid[0] ? 1.0 : -1.0;
  } else {
    // Adapted from https://mfem.org/howto/outer_normals/
    mfem::FaceElementTransformations face_element_transformations;
    mfem::IsoparametricTransformation inner_element_transformation, outer_element_transformation;
    mesh.GetFaceElementTransformations(
      element_face_to_mesh_face[element_face],
      face_element_transformations,
      inner_element_transformation,
      outer_element_transformation
    );
    face_element_transformations.SetIntPoint(&mfem::Geometries.GetCenter(face_element_transformations.GetGeometryType()));
    mfem::CalcOrtho(face_element_transformations.Jacobian(), unit_normal);
    unit_normal /= unit_normal.Norml2();

    if (face_element_transformations.Elem2No == element) {
      unit_normal *= -1.0;
    }
  }

  return unit_normal;
}

std::tuple<int, int, bool> getElementFaceOfBoundaryElement(const mfem::Mesh& mesh, int boundary_element) {
  assert(boundary_element >= 0);
  assert(boundary_element < mesh.GetNBE());

  int element = -1;
  int element_face = -1;
  bool element_face_found = false;
  if (mesh.GetBdrElementFaceIndex(boundary_element) < mesh.GetNumFaces()) {
    int element_info;
    mesh.GetBdrElementAdjacentElement(boundary_element, element, element_info);
    mfem::Array<int> element_face_to_mesh_face = elementFaceToMeshFaceMap(mesh, element);
    for (element_face = 0; element_face < element_face_to_mesh_face.Size(); element_face++) {
      if (mesh.GetBdrElementFaceIndex(boundary_element) == element_face_to_mesh_face[element_face]) {
        element_face_found = true;
        break;
      }
    }
  }

  return std::make_tuple(element, element_face, element_face_found);
}

double getMeshVolume(mfem::Mesh& mesh) {
  double mesh_volume = 0.0;
  for (int element = 0; element < mesh.GetNE(); element++) {
    mesh_volume += mesh.GetElementVolume(element);
  }
  return mesh_volume;
}

int getNumFacesOnElement(const mfem::Mesh& mesh, int element) {
  assert(element < mesh.GetNE());
  const mfem::Geometry::Type element_geometry_type = mesh.GetElementGeometry(element);
  const mfem::Geometry geometry;
  return geometry.NumBdr(element_geometry_type);
}

std::unordered_map<std::string, int> getSideNameToBoundaryAttributeForInlineMeshes(const int mesh_dimension) {
  assert(mesh_dimension > 0);
  assert(mesh_dimension < 4);

  std::unordered_map<std::string, int> side_name_to_boundary_attribute;

  if (mesh_dimension == 1) {
    side_name_to_boundary_attribute["left"] = 1;
    side_name_to_boundary_attribute["right"] = 2;
  } else if (mesh_dimension == 2) {
    side_name_to_boundary_attribute["bottom"] = 1;
    side_name_to_boundary_attribute["right"] = 2;
    side_name_to_boundary_attribute["top"] = 3;
    side_name_to_boundary_attribute["left"] = 4;
  } else if (mesh_dimension == 3) {
    side_name_to_boundary_attribute["back"] = 1;
    side_name_to_boundary_attribute["bottom"] = 2;
    side_name_to_boundary_attribute["right"] = 3;
    side_name_to_boundary_attribute["top"] = 4;
    side_name_to_boundary_attribute["left"] = 5;
    side_name_to_boundary_attribute["front"] = 6;
  }

  return side_name_to_boundary_attribute;
}

double getSmallestCellLengthscale(mfem::Mesh& mesh) {
  double min_length = std::numeric_limits<double>::max();
  for (int i = 0; i < mesh.GetNE(); ++i)
    min_length = fmin(min_length, mesh.GetElementSize(i, 1));

  return min_length;
}

} // namespace mfpic
