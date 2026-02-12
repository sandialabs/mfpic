#include <libmfpic/Errors.hpp>
#include <libmfpic/MeshUtilities.hpp>
#include <libmfpic/Particle.hpp>
#include <libmfpic/PeriodicParticleBoundary.hpp>

namespace mfpic {

ElementFaceContainer<std::shared_ptr<ParticleBoundary>> PeriodicParticleBoundary::generatePeriodicParticleBoundaries(mfem::Mesh& mesh) {
  auto periodic_boundary = std::make_shared<PeriodicParticleBoundary>();
  ElementFaceContainer<std::shared_ptr<ParticleBoundary>> periodic_boundaries;

  constexpr double periodic_face_consistency_check_relative_tolerance = 1.0e-2;
  const int dimension = mesh.Dimension();
  const double characteristic_element_size = std::pow(mesh.GetElementVolume(0), 1./dimension);
  const double periodic_face_consistency_check_absolute_tolerance = periodic_face_consistency_check_relative_tolerance * characteristic_element_size;
  for (int face = 0; face < mesh.GetNumFaces(); face++) {
    const bool face_is_periodic = mesh.GetFaceElementTransformations(face)->CheckConsistency() > periodic_face_consistency_check_absolute_tolerance;
    if (face_is_periodic) {
      mfem::FaceElementTransformations face_element_transformations;
      mfem::IsoparametricTransformation& element_1_transformation = periodic_boundary->element_transformation_dump_.emplace_front();
      mfem::IsoparametricTransformation& element_2_transformation = periodic_boundary->element_transformation_dump_.emplace_front();
      mesh.GetFaceElementTransformations(
        face,
        face_element_transformations,
        element_1_transformation,
        element_2_transformation
      );
      const int elements_attached_to_this_face[] = {face_element_transformations.Elem1No, face_element_transformations.Elem2No};
      for (int i = 0; i < 2; i++) {
        const int element = elements_attached_to_this_face[i];
        const int other_element = elements_attached_to_this_face[1-i];
        for (int local_face = 0; local_face < getNumFacesOnElement(mesh, element); local_face++) {
          if (getElementOnOtherSideOfFace(mesh, element, local_face) == other_element) {
            periodic_boundary->face_element_transformations_.insert(element, local_face, face_element_transformations);
            periodic_boundaries.insert(element, local_face, periodic_boundary);
            break;
          }
        }
      }
    }
  }

  return periodic_boundaries;
}

Particle PeriodicParticleBoundary::applyBoundary(int element_face, const Particle& incoming_particle) const {
  Particle translated_particle = incoming_particle;

  const int current_element = incoming_particle.element;
  mfem::FaceElementTransformations face_element_transformations = face_element_transformations_.at(
    current_element,
    element_face
  );
  const bool current_element_is_face_element_1 = current_element == face_element_transformations.Elem1No;
  mfem::IsoparametricTransformation current_element_transformation = dynamic_cast<mfem::IsoparametricTransformation&>(
    current_element_is_face_element_1 ? face_element_transformations.GetElement1Transformation() : face_element_transformations.GetElement2Transformation());
  mfem::IntegrationPointTransformation face_to_current_element_transformation =
    current_element_is_face_element_1 ? face_element_transformations.GetIntPoint1Transformation() : face_element_transformations.GetIntPoint2Transformation();
  mfem::IsoparametricTransformation new_element_transformation = dynamic_cast<mfem::IsoparametricTransformation&>(
    current_element_is_face_element_1 ? face_element_transformations.GetElement2Transformation() : face_element_transformations.GetElement1Transformation());
  mfem::IntegrationPointTransformation face_to_new_element_transformation =
    current_element_is_face_element_1 ? face_element_transformations.GetIntPoint2Transformation() : face_element_transformations.GetIntPoint1Transformation();

  mfem::IntegrationPoint face_reference_coords;
  if (face_element_transformations.GetGeometryType() == mfem::Geometry::POINT) {
    face_reference_coords.x = 0.0;
  } else {
    const int mesh_dimensions = face_element_transformations.GetSpaceDim();

    mfem::IntegrationPoint current_element_reference_coords;
    mfem::Vector incoming_particle_position(incoming_particle.position.GetData(), mesh_dimensions);
    current_element_transformation.TransformBack(incoming_particle_position, current_element_reference_coords);
    mfem::Vector current_element_reference_coords_as_vector(mesh_dimensions);
    switch (mesh_dimensions) {
    case 3:
      current_element_reference_coords_as_vector[2] = current_element_reference_coords.z;
      [[fallthrough]];
    case 2:
      current_element_reference_coords_as_vector[1] = current_element_reference_coords.y;
      [[fallthrough]];
    case 1:
      current_element_reference_coords_as_vector[0] = current_element_reference_coords.x;
      break;
    default:
      errorWithDeveloperMessage("Invalid mesh dimensions.");
    }

    face_to_current_element_transformation.Transf.TransformBack(current_element_reference_coords_as_vector, face_reference_coords);
  }

  mfem::IntegrationPoint new_element_reference_coords;
  face_to_new_element_transformation.Transform(face_reference_coords, new_element_reference_coords);

  mfem::Vector new_element_physical_coords;
  new_element_transformation.Transform(new_element_reference_coords, new_element_physical_coords);

  for (int i = 0; i < new_element_physical_coords.Size(); i++) {
    translated_particle.position[i] = new_element_physical_coords[i];
  }

  translated_particle.element =
    current_element_is_face_element_1 ? face_element_transformations.Elem2No : face_element_transformations.Elem1No;

  return translated_particle;
}

} // namespace mfpic
