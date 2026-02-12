#include <libmfpic/IntegratedCharge.hpp>
#include <libmfpic/MeshUtilities.hpp>
#include <libmfpic/ParticleContainer.hpp>
#include <libmfpic/ParticleOperations.hpp>
#include <libmfpic/PeriodicParticleBoundary.hpp>

#include <mfem/mfem.hpp>

#include <limits>

namespace mfpic {

ParticleOperations::ParticleOperations(
  Discretization &discretization,
  std::vector<std::shared_ptr<ParticleBoundaryFactory>> particle_boundary_factories,
  std::shared_ptr<ParticleBoundaryFactory> default_particle_boundary_factory
) :
  discretization_(discretization),
  dim_(discretization_.getFeSpace().GetMesh()->Dimension())
{
  mfem::Mesh& mesh = *discretization_.getFeSpace().GetMesh();
  element_face_unit_normal_ = std::make_shared<ElementFaceContainer<mfem::Vector>>();
  for (int element = 0; element < mesh.GetNE(); element++) {
    const int num_faces = getNumFacesOnElement(mesh, element);
    num_faces_on_element_.push_back(num_faces);

    for (int face = 0; face < num_faces; face++) {
      mfem::Vector face_unit_normal = getElementFaceOutwardUnitNormal(mesh, element, face);
      element_face_unit_normal_->insert(element, face, face_unit_normal);
      mfem::Vector face_centroid = getElementFaceCentroid(mesh, element, face);
      element_face_centroid_dot_unit_normal_.insert(element, face, face_centroid * face_unit_normal);
      element_face_other_element_.insert(element, face, getElementOnOtherSideOfFace(mesh, element, face));
    }
  }

  ParticleBoundaryFactory::Parameters particle_boundary_factory_params{element_face_unit_normal_};

  std::shared_ptr<ParticleBoundary> default_particle_boundary = default_particle_boundary_factory->createBoundary(
    particle_boundary_factory_params
  );

  std::unordered_map<int, std::shared_ptr<ParticleBoundary>> attribute_to_boundary;
  for (std::shared_ptr<ParticleBoundaryFactory> particle_boundary_factory : particle_boundary_factories) {
    const auto emplace_result = attribute_to_boundary.emplace(
      particle_boundary_factory->getBoundaryAttribute(),
      particle_boundary_factory->createBoundary(particle_boundary_factory_params)
    );
    [[maybe_unused]] const bool emplace_succeeded = emplace_result.second;
    assert(emplace_succeeded);
  }

  particle_boundaries_ = PeriodicParticleBoundary::generatePeriodicParticleBoundaries(mesh);
  for (int boundary_element = 0; boundary_element < mesh.GetNBE(); boundary_element++) {
    const int attribute = mesh.GetBdrAttribute(boundary_element);
    const auto [element, element_face, element_face_exists] = getElementFaceOfBoundaryElement(mesh, boundary_element);
    if (not element_face_exists) continue;
    if (attribute_to_boundary.contains(attribute)) {
      particle_boundaries_.insert(element, element_face, attribute_to_boundary.at(attribute));
    } else {
      particle_boundaries_.insert(element, element_face, default_particle_boundary);
    }
  }
}

ParticleContainer ParticleOperations::accelerate(
  double dt,
  const ParticleContainer& current_particles,
  const ElectromagneticFieldsEvaluator& field_provider
) const {
  ParticleContainer accelerated_particles = current_particles;

  for (Particle& particle : accelerated_particles) {
    if (not particle.is_alive) continue;

    const mfem::Vector position(particle.position.GetData(), dim_);

    particle.velocity.Add(
      dt * particle.species.charge_over_mass,
      field_provider.getEFieldAt(position, particle.element)
    );
  }

  return accelerated_particles;
}

ParticleContainer ParticleOperations::move(
  double dt,
  const ParticleContainer& current_particles
) const {
  ParticleContainer moved_particles = current_particles;

  const ElementFaceContainer<mfem::Vector>& element_face_unit_normal = *element_face_unit_normal_;
  for (Particle& particle : moved_particles) {
    if (not particle.is_alive) continue;

    double time_remaining = dt;
    do {
      const int current_element = particle.element;
      const int num_faces = num_faces_on_element_[current_element];

      const mfem::Vector position(particle.position.GetData(), dim_);
      const mfem::Vector velocity(particle.velocity.GetData(), dim_);

      int closest_face = -1;
      double time_to_closest_face = std::numeric_limits<double>::max();
      for (int face = 0; face < num_faces; face++) {
        const double distance_to_face =
          element_face_centroid_dot_unit_normal_.at(current_element, face) -
          position * element_face_unit_normal.at(current_element, face);
        const double speed_toward_face = velocity * element_face_unit_normal.at(current_element, face);
        const double time_to_face = speed_toward_face > 0.0 ? distance_to_face / speed_toward_face : std::numeric_limits<double>::max();
        const bool particle_will_cross_this_face_first = 0.0 <= time_to_face and time_to_face < time_to_closest_face;
        if (particle_will_cross_this_face_first) {
          closest_face = face;
          time_to_closest_face = time_to_face;
        }
      }

      const double time_spent_in_element = std::min(time_to_closest_face, time_remaining);
      const bool particle_is_crossing_a_face = time_to_closest_face <= time_remaining;
      particle.position.Add(time_spent_in_element, particle.velocity);
      time_remaining -= time_spent_in_element;
      if (particle_is_crossing_a_face) {
        if (particle_boundaries_.contains(current_element, closest_face)) {
          particle = particle_boundaries_.at(current_element, closest_face)->applyBoundary(closest_face, particle);
          if (not particle.is_alive) break;
        } else {
          particle.element = element_face_other_element_.at(current_element, closest_face);
          assert(particle.element >= 0);
        }
      }
    } while (time_remaining > 0.0);
  }

  return moved_particles;
}

IntegratedCharge ParticleOperations::assembleCharge(
  const ParticleContainer& current_particles
) const {
  ParticleContainer particles = current_particles;
  mfem::IntegrationPoint integration_point;
  mfem::Array<int> vector_dofs;
  mfem::FiniteElementSpace finite_element_space = discretization_.getFeSpace();

  IntegratedCharge charge_state(discretization_);

  charge_state.setIntegratedChargeValue(0.0);
  mfem::Mesh &mesh = *finite_element_space.GetMesh();

  for (Particle& particle : particles) {
    if (not particle.is_alive) continue;

    const int elem_id = particle.element;
    const double particle_charge = particle.species.charge;
    mfem::ElementTransformation * element_transformation = mesh.GetElementTransformation(elem_id);
    const mfem::FiniteElement *fe = finite_element_space.GetFE(elem_id);

    const mfem::Vector particle_position(particle.position.GetData(), dim_); 
    element_transformation->TransformBack(particle_position, integration_point);
    element_transformation->SetIntPoint(&integration_point);
    mfem::Vector psi_i(fe->GetDof());
    fe->CalcPhysShape(*element_transformation,psi_i);
    finite_element_space.GetElementVDofs(elem_id, vector_dofs);

    for (int i = 0; i < fe->GetDof(); i++) {
      charge_state.addIntegratedChargeValue(vector_dofs[i],particle.weight * particle_charge * psi_i(i));
    }
  }

  return charge_state;
  }

} // namespace mfpic
