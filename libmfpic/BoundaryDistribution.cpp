#include <libmfpic/BoundaryDistribution.hpp>
#include <libmfpic/Errors.hpp>
#include <libmfpic/MeshUtilities.hpp>

namespace mfpic {

namespace {

double computeBoundaryElementArea(mfem::Mesh& mesh, int boundary_element) {
  if (mesh.Dimension() == 1) return 1.0;

  mfem::ElementTransformation* transformation = mesh.GetBdrElementTransformation(boundary_element);
  const mfem::IntegrationRule& rule = mfem::IntRules.Get(
    mesh.GetBdrElementGeometry(boundary_element),
    2 * transformation->OrderJ());
  double area = 0.0;
  for (int i_point = 0; i_point < rule.GetNPoints(); i_point++) {
    const mfem::IntegrationPoint& point = rule.IntPoint(i_point);
    transformation->SetIntPoint(&point);
    area += point.weight * transformation->Weight();
  }
  return area;
}

mfem::Vector padToLength3(const mfem::Vector& vector) {
  mfem::Vector padded({0.0, 0.0, 0.0});
  for (int i = 0; i < vector.Size(); i++) {
    padded[i] = vector[i];
  }
  return padded;
}

} // namespace

BoundaryDistribution::BoundaryDistribution(std::shared_ptr<mfem::Mesh> mesh, int boundary_attribute) :
  mesh_(mesh)
{
  std::vector<double> face_areas;
  for (int boundary_element = 0; boundary_element < mesh_->GetNBE(); boundary_element++) {
    if (mesh_->GetBdrAttribute(boundary_element) != boundary_attribute) continue;
    const auto [element, element_face, element_face_exists] = getElementFaceOfBoundaryElement(*mesh_, boundary_element);
    if (not element_face_exists) continue;

    mfem::Vector outward_unit_normal = padToLength3(getElementFaceOutwardUnitNormal(*mesh_, element, element_face));
    faces_info_.push_back(FaceInfo_{
      .boundary_element = boundary_element,
      .element = element,
      .element_face = element_face,
      .outward_unit_normal = outward_unit_normal,
    });
    face_areas.push_back(computeBoundaryElementArea(*mesh_, boundary_element));
  }

  if (faces_info_.empty()) {
    errorWithUserMessage(
      "Mesh side with boundary attribute " + std::to_string(boundary_attribute) +
      " has no non-periodic boundary faces!");
  }

  total_area_ = std::reduce(face_areas.begin(), face_areas.end());
  face_distribution_ = std::discrete_distribution<>(face_areas.begin(), face_areas.end());
}

double BoundaryDistribution::getTotalArea() const {
  return total_area_;
}

BoundaryPoint BoundaryDistribution::generateRandomBoundaryPoint(RandomNumberGenerator& generator) {
  const FaceInfo_& random_face_info = faces_info_[face_distribution_(generator)];

  mfem::Vector physical_point;
  if (mesh_->Dimension() == 1) {
    physical_point = getElementFaceCentroid(*mesh_, random_face_info.element, random_face_info.element_face);
  } else {
    mfem::IntegrationPoint random_reference_point;
    // TODO MFEM currently calls rand() to generate random numbers. oughtta use c++ generators.
    // The following is done to provide a seed consistent with our own.
    std::uniform_int_distribution<unsigned> random_unsigned_distribution;
    std::srand(random_unsigned_distribution(generator));
    mfem::Geometry::GetRandomPoint(mesh_->GetBdrElementGeometry(random_face_info.boundary_element), random_reference_point);

    mfem::ElementTransformation* transformation = mesh_->GetBdrElementTransformation(random_face_info.boundary_element);
    transformation->Transform(random_reference_point, physical_point);
  }

  return BoundaryPoint{
    .position = padToLength3(physical_point),
    .element = random_face_info.element,
    .outward_unit_normal = random_face_info.outward_unit_normal,
  };
}

} // namespace mfpic
