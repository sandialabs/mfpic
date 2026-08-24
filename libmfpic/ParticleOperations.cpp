#include <libmfpic/BuildVarianceReductionParametersFromYaml.hpp>
#include <libmfpic/DGEulerOperations.hpp>
#include <libmfpic/Constants.hpp>
#include <libmfpic/IntegratedCharge.hpp>
#include <libmfpic/LowFidelityOperations.hpp>
#include <libmfpic/LowFidelityState.hpp>
#include <libmfpic/MeshUtilities.hpp>
#include <libmfpic/ParticleContainer.hpp>
#include <libmfpic/ParticleOperations.hpp>
#include <libmfpic/PeriodicParticleBoundary.hpp>
#include <libmfpic/Species.hpp>

#include <mfem/linalg/vector.hpp>
#include <mfem/mfem.hpp>

#include <limits>
#include <unordered_map>

namespace mfpic {

ParticleOperations::ParticleOperations(
  Discretization &discretization,
  std::vector<std::shared_ptr<ParticleBoundaryFactory>> particle_boundary_factories,
  std::shared_ptr<ParticleBoundaryFactory> default_particle_boundary_factory,
  std::unordered_map<std::string, Species> species_map,
  const int velocity_dims
) :
  discretization_(discretization),
  dim_(discretization_.getFeSpace().GetMesh()->Dimension()),
  velocity_dims_(velocity_dims)
{
  mfem::Mesh& mesh = *discretization_.getFeSpace().GetMesh();

  for (const auto & [name, species] : species_map) {
    particle_number_density_.insert({species, mfem::Vector(mesh.GetNE())});
    particle_bulk_velocity_.insert({species, mfem::DenseMatrix(3, mesh.GetNE())});
    particle_temperature_.insert({species, mfem::Vector(mesh.GetNE())});
    variance_reduced_particle_number_density_.insert({species, mfem::Vector(mesh.GetNE())});
    variance_reduced_particle_bulk_velocity_.insert({species, mfem::DenseMatrix(3, mesh.GetNE())});
    variance_reduced_particle_temperature_.insert({species, mfem::Vector(mesh.GetNE())});
    avg_particle_number_density_.insert({species, 0.0});
    avg_particle_bulk_velocity_.insert({species, mfem::Vector(3)});
    avg_particle_temperature_.insert({species, 0.0});
    avg_variance_reduced_particle_number_density_.insert({species, 0.0});
    avg_variance_reduced_particle_bulk_velocity_.insert({species, mfem::Vector(3)});
    avg_variance_reduced_particle_temperature_.insert({species, 0.0});
    variance_reduced_postprocessors_.noise_reducing_factor.insert({species,mfem::Vector(mesh.GetNE())});
    variance_reduced_postprocessors_.particle_distribution_function.insert({species,mfem::Vector(mesh.GetNE())});
    variance_reduced_postprocessors_.low_fidelity_particle_distribution_function.insert({species,mfem::Vector(mesh.GetNE())});
    sum_of_weights_.insert({species, mfem::Vector(mesh.GetNE())});
  }

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

  #pragma omp parallel for
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
  #pragma omp parallel for
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

  for (const Particle& particle : particles) {
    if (not particle.is_alive) continue;

    const int elem_id = particle.element;
    const double particle_charge = particle.species.charge;
    mfem::ElementTransformation * element_transformation = mesh.GetElementTransformation(elem_id);
    const mfem::FiniteElement *fe = finite_element_space.GetFE(elem_id);

    const mfem::Vector particle_position(particle.position.GetData(), dim_); 
    element_transformation->TransformBack(particle_position, integration_point);
    element_transformation->SetIntPoint(&integration_point);
    mfem::Vector psi(fe->GetDof());
    fe->CalcPhysShape(*element_transformation,psi);
    finite_element_space.GetElementVDofs(elem_id, vector_dofs);

    for (int i = 0; i < fe->GetDof(); i++) {
      charge_state.addIntegratedChargeValue(vector_dofs[i],particle.weight * particle_charge * psi(i));
    }
  }

  return charge_state;
}

IntegratedCharge ParticleOperations::assembleVarianceReducedCharge(
  const ParticleContainer& current_particles,
  const LowFidelityState& low_fidelity_state,
  const LowFidelityOperations& low_fidelity_operations
) const {
  ParticleContainer particles = current_particles;

  mfem::IntegrationPoint integration_point;
  mfem::Array<int> vector_dofs;

  mfem::FiniteElementSpace finite_element_space = discretization_.getFeSpace();
  mfem::Mesh& mesh = *finite_element_space.GetMesh();

  std::unordered_map<Species, std::vector<char>> variance_reduction_performed; 
  variance_reduction_performed.reserve(variance_reduced_particle_number_density_.size());
  for (const auto& kv : variance_reduced_particle_number_density_)
    variance_reduction_performed.emplace(kv.first, std::vector<char>(finite_element_space.GetNDofs(), 0));

  mfem::Vector particle_velocity(velocity_dims_);  
  for (const Particle& particle : particles) {
    if (!particle.is_alive) continue;

    for (int vel_dim = 0; vel_dim < velocity_dims_; ++vel_dim)
      particle_velocity(vel_dim) = particle.velocity(vel_dim);

    const int element_id = particle.element;
    const mfem::FiniteElement* finite_element = finite_element_space.GetFE(element_id);
    finite_element_space.GetElementVDofs(element_id, vector_dofs);

    const mfem::Vector particle_position(particle.position.GetData(), dim_);
    const int low_fidelity_species_index = low_fidelity_state.getSpeciesIndex(particle.species);
    if (low_fidelity_species_index < 0) continue;

    const double low_fidelity_value =
      low_fidelity_operations.evaluateParticleDistributionFunction(
        low_fidelity_state, particle_position, particle_velocity, particle.element, low_fidelity_species_index);

    const double noise_reducing_factor =
      1.0 - low_fidelity_value / particle.particle_distribution_function_value;

    if (std::abs(noise_reducing_factor) < 1.0) {
      for (int local_dof = 0; local_dof < finite_element->GetDof(); ++local_dof) {
        const int global_dof = vector_dofs[local_dof];
        variance_reduction_performed.at(particle.species)[global_dof] = 1;
      }
    }
  }

  IntegratedCharge integrated_charge(discretization_);
  integrated_charge.setIntegratedChargeValue(0.0);

  for (const Particle& particle : particles) {
    if (!particle.is_alive) continue;

    for (int vel_dim = 0; vel_dim < velocity_dims_; ++vel_dim)
      particle_velocity(vel_dim) = particle.velocity(vel_dim);

    const int element_id = particle.element;
    const Species& particle_species = particle.species;
    const double particle_charge = particle_species.charge;

    mfem::ElementTransformation* element_transformation =
      mesh.GetElementTransformation(element_id);

    const mfem::FiniteElement* finite_element =
      finite_element_space.GetFE(element_id);

    const mfem::Vector particle_position(particle.position.GetData(), dim_);

    element_transformation->TransformBack(particle_position, integration_point);
    element_transformation->SetIntPoint(&integration_point);

    mfem::Vector shape_functions(finite_element->GetDof());
    finite_element->CalcPhysShape(*element_transformation, shape_functions);

    finite_element_space.GetElementVDofs(element_id, vector_dofs);
    for (int local_dof = 0; local_dof < finite_element->GetDof(); ++local_dof) {
      const int global_dof = vector_dofs[local_dof];
      double noise_reducing_factor = 1.0;
      if (variance_reduction_performed.at(particle.species)[global_dof] == 1)
      {
        const int low_fidelity_species_index = low_fidelity_state.getSpeciesIndex(particle.species);
        const double low_fidelity_value = low_fidelity_operations.evaluateParticleDistributionFunction(low_fidelity_state, particle_position, particle_velocity, element_id, low_fidelity_species_index);
        noise_reducing_factor = 1.0 - low_fidelity_value / particle.particle_distribution_function_value;
      }
      integrated_charge.addIntegratedChargeValue(global_dof, particle.weight * particle_charge * shape_functions(local_dof) * noise_reducing_factor);
    }
  }

  for (int ispecies = 0; ispecies < low_fidelity_state.numSpecies(); ++ispecies) {
    const Species species = low_fidelity_state.getSpeciesState(ispecies).getSpecies();
    IntegratedCharge low_fidelity_charge_state = low_fidelity_operations.assembleChargePerSpecies(low_fidelity_state,ispecies);
    for (int global_dof = 0; global_dof < finite_element_space.GetNDofs(); ++global_dof) {
      if (variance_reduction_performed.at(species)[global_dof] == 1) {
        integrated_charge.addIntegratedChargeValue(
        global_dof,
        low_fidelity_charge_state.getIntegratedChargeValue(global_dof));
      }
    }
  }
  return integrated_charge;
}

std::unordered_map<Species, mfem::Vector>& ParticleOperations::getNumberDensity(const ParticleContainer& particles
) {

  for (auto & species_and_number_density : particle_number_density_)
    species_and_number_density.second = 0.0;

  mfem::Array<int> vector_dofs;
  mfem::FiniteElementSpace finite_element_space = discretization_.getFeSpace();
  mfem::Mesh &mesh = *finite_element_space.GetMesh();

  for (const Particle& particle : particles) {
    if (not particle.is_alive) continue;

    const int elem_id = particle.element;
    const Species & species = particle.species;
    particle_number_density_.at(species)(elem_id) += particle.weight / mesh.GetElementVolume(elem_id);
  }

  return this->particle_number_density_;
}

std::unordered_map<Species, double>& ParticleOperations::getAvgNumberDensity(const ParticleContainer& particles
) {

  for (auto & species_and_number_density : avg_particle_number_density_)
    species_and_number_density.second = 0.0;

  mfem::Array<int> vector_dofs;
  mfem::FiniteElementSpace finite_element_space = discretization_.getFeSpace();
  mfem::Mesh &mesh = *finite_element_space.GetMesh();

  for (const Particle& particle : particles) {
    if (not particle.is_alive) continue;

    const Species & species = particle.species;
    avg_particle_number_density_.at(species) += particle.weight / getMeshVolume(mesh);
  }

  return this->avg_particle_number_density_;
}


std::unordered_map<Species, mfem::Vector>& ParticleOperations::getVarianceReducedNumberDensity(
  const ParticleContainer& particles, 
  const LowFidelityState& low_fidelity_state,
  const DGEulerOperations& low_fidelity_operations
) {

  for (auto & species_and_number_density : variance_reduced_particle_number_density_)
    species_and_number_density.second = 0.0;

  mfem::FiniteElementSpace finite_element_space = discretization_.getFeSpace();
  mfem::Mesh &mesh = *finite_element_space.GetMesh();

  std::unordered_map<Species,mfem::Vector> low_fidelity_integral = low_fidelity_operations.integralForVarianceReducedNumberDensity(finite_element_space, low_fidelity_state);
  this->getVarianceReducedPostprocessors(particles,low_fidelity_state,low_fidelity_operations);

  std::unordered_map<Species, std::vector<char>> variance_reduction_performed; 
  variance_reduction_performed.reserve(variance_reduced_particle_number_density_.size());
  for (const auto& kv : variance_reduced_particle_number_density_)
    variance_reduction_performed.emplace(kv.first, std::vector<char>(finite_element_space.GetNE(), 0));

  /////
  //Compute standard PIC moments
  std::unordered_map<Species, mfem::Vector> number_density     = this->getNumberDensity(particles);
  std::unordered_map<Species, mfem::DenseMatrix> bulk_velocity = this->getBulkVelocity(particles,true);
  std::unordered_map<Species, mfem::Vector> temperature        = this->getTemperature(particles,true,true);
  /////

  mfem::Vector particle_velocity(velocity_dims_);  
  for (const Particle& particle : particles) {
    if (not particle.is_alive) continue;

    for (int vel_dim = 0; vel_dim < velocity_dims_; ++vel_dim)
      particle_velocity(vel_dim) = particle.velocity(vel_dim);

    const int elem_id = particle.element;
    const double element_volume = mesh.GetElementVolume(elem_id);

    const mfem::Vector particle_position(particle.position.GetData(), dim_);

    const int low_fidelity_species_index =
      low_fidelity_state.getSpeciesIndex(particle.species);
    
    bool perform_variance_reduction = (low_fidelity_species_index >= 0);
    if (variance_reduction_parameters_.limit_variance_reduction)
    {
      perform_variance_reduction = (low_fidelity_species_index >= 0) &&
        (variance_reduced_postprocessors_.noise_reducing_factor
            .at(particle.species)(elem_id) < 1.0);
    }

    if (perform_variance_reduction)
    {
      variance_reduction_performed.at(particle.species)[elem_id] = 1;

      double low_fidelity_particle_distribution_function_value = 0.0;
      if (variance_reduction_parameters_.strategy == VarianceReductionParameters::Strategy::LocalMaxwellian)
      {
        mfem::Vector primitive_state(5);
        primitive_state(0) = number_density.at(particle.species)[elem_id];
        const mfem::Vector bulk_velocity_in_element(bulk_velocity.at(particle.species).GetColumn(elem_id), 3);
        primitive_state(0) = number_density.at(particle.species)(elem_id);
        primitive_state(1) = bulk_velocity_in_element(0);
        primitive_state(2) = bulk_velocity_in_element(1);
        primitive_state(3) = bulk_velocity_in_element(2);
        primitive_state(4) = temperature.at(particle.species)(elem_id);
        low_fidelity_particle_distribution_function_value = euler::evaluateMaxwellian(primitive_state,particle_velocity,particle.species);
      }
      else
        low_fidelity_particle_distribution_function_value = low_fidelity_operations.evaluateParticleDistributionFunction(low_fidelity_state,particle_position,particle_velocity,particle.element,low_fidelity_species_index);

      double noise_reducing_factor = (1 - low_fidelity_particle_distribution_function_value / particle.particle_distribution_function_value);
      // std::cout << "noise reducing factor: " << noise_reducing_factor << std::endl;
      // std::cout << "fnr : " << low_fidelity_particle_distribution_function_value << std::endl;
      // std::cout << "f: " << particle.particle_distribution_function_value << std::endl;
      variance_reduced_particle_number_density_.at(particle.species)(elem_id) += (particle.weight * noise_reducing_factor) / element_volume;
    }
    else
    {
      variance_reduced_particle_number_density_.at(particle.species)(elem_id) += particle.weight / element_volume;
    }
  }
  for (int elem_id = 0; elem_id < finite_element_space.GetNE(); ++elem_id) {
    const double element_volume = mesh.GetElementVolume(elem_id);
    for (int ispecies = 0; ispecies < low_fidelity_state.numSpecies(); ++ispecies) {
      const Species species = low_fidelity_state.getSpeciesState(ispecies).getSpecies();
      if (variance_reduction_performed.at(species)[elem_id] == 1) {
        if (variance_reduction_parameters_.specified_lf_moments)
        {
          variance_reduced_particle_number_density_.at(species)(elem_id) += variance_reduction_parameters_.reference_number_density;
        }
        else
        {
          variance_reduced_particle_number_density_.at(species)(elem_id) +=
            low_fidelity_integral.at(species)(elem_id) / element_volume;
        }
      }
    }
  }
  return this->variance_reduced_particle_number_density_;
}

std::unordered_map<Species, mfem::DenseMatrix>& ParticleOperations::getBulkVelocity(const ParticleContainer& particles, const bool sum_weights
) {

  for (auto & species_and_bulk_velocity : particle_bulk_velocity_)
    species_and_bulk_velocity.second = 0.0;

  if (sum_weights) this->sumParticleWeights_(particles);

  for (const Particle& particle : particles) {
    if (not particle.is_alive) continue;

    const int elem_id = particle.element;
    const Species & species = particle.species;
    const double sum_weights = sum_of_weights_.at(species)(elem_id);
    if (sum_weights <= 0.0) continue;

    mfem::Vector velocity_in_element(particle_bulk_velocity_.at(species).GetColumn(elem_id), 3);
    velocity_in_element.Add(particle.weight / sum_weights, particle.velocity);
  }

  return this->particle_bulk_velocity_;
}

std::unordered_map<Species, mfem::Vector>& ParticleOperations::getAvgBulkVelocity(const ParticleContainer& particles, const bool sum_weights
) {

  for (auto & species_and_bulk_velocity : avg_particle_bulk_velocity_)
    species_and_bulk_velocity.second = 0.0;

  if (sum_weights) this->sumParticleWeights_(particles);

  for (const Particle& particle : particles) {
    if (not particle.is_alive) continue;

    const int elem_id = particle.element;
    const Species & species = particle.species;
    const double sum_weights = sum_of_weights_.at(species)(elem_id);
    const double number_of_physical_particles = sum_of_weights_.at(particle.species).Sum();
    if (sum_weights <= 0.0) continue;

    avg_particle_bulk_velocity_.at(particle.species)(0) = particle.weight * particle.velocity(0) / number_of_physical_particles;
    avg_particle_bulk_velocity_.at(particle.species)(1) = particle.weight * particle.velocity(1) / number_of_physical_particles;
    avg_particle_bulk_velocity_.at(particle.species)(2) = particle.weight * particle.velocity(2) / number_of_physical_particles;
  }

  return this->avg_particle_bulk_velocity_;
}

std::unordered_map<Species, mfem::DenseMatrix>& ParticleOperations::getVarianceReducedBulkVelocity(
  const ParticleContainer& particles, 
  const LowFidelityState& low_fidelity_state,
  const DGEulerOperations& low_fidelity_operations
) {

  for (auto & species_and_bulk_velocity : variance_reduced_particle_bulk_velocity_)
    species_and_bulk_velocity.second = 0.0;

  this->sumParticleWeights_(particles);
  this->getVarianceReducedPostprocessors(particles,low_fidelity_state,low_fidelity_operations);
  mfem::FiniteElementSpace finite_element_space = discretization_.getFeSpace();
  mfem::Mesh &mesh = *finite_element_space.GetMesh();
  std::unordered_map<Species, mfem::DenseMatrix> low_fidelity_integral = low_fidelity_operations.integralForVarianceReducedBulkVelocity(finite_element_space, low_fidelity_state);

  std::unordered_map<Species, std::vector<char>> variance_reduction_performed; 
  variance_reduction_performed.reserve(variance_reduced_particle_number_density_.size());
  for (const auto& kv : variance_reduced_particle_number_density_)
    variance_reduction_performed.emplace(kv.first, std::vector<char>(finite_element_space.GetNE(), 0));

  /////
  //Compute standard PIC moments
  std::unordered_map<Species, mfem::Vector> number_density     = this->getNumberDensity(particles);
  std::unordered_map<Species, mfem::DenseMatrix> bulk_velocity = this->getBulkVelocity(particles,true);
  std::unordered_map<Species, mfem::Vector> temperature        = this->getTemperature(particles,true,true);
  /////

  mfem::Vector particle_velocity(velocity_dims_);  
  for (const Particle& particle : particles) {
    if (not particle.is_alive) continue;

    for (int vel_dim = 0; vel_dim < velocity_dims_; ++vel_dim)
      particle_velocity(vel_dim) = particle.velocity(vel_dim);

    const int elem_id = particle.element;
    const double element_volume = mesh.GetElementVolume(elem_id);

    const mfem::Vector particle_position(particle.position.GetData(), dim_);
    const int low_fidelity_species_index =
      low_fidelity_state.getSpeciesIndex(particle.species);

    bool perform_variance_reduction = (low_fidelity_species_index >= 0);
    if (variance_reduction_parameters_.limit_variance_reduction)
    {
      perform_variance_reduction = (low_fidelity_species_index >= 0) &&
        (variance_reduced_postprocessors_.noise_reducing_factor
            .at(particle.species)(elem_id) < 1.0);
    }

    mfem::Vector velocity_in_element(variance_reduced_particle_bulk_velocity_.at(particle.species).GetColumn(elem_id), 3);
    double variance_reduced_number_density = variance_reduced_particle_number_density_.at(particle.species)(elem_id);
    if (perform_variance_reduction)
    {

      variance_reduction_performed.at(particle.species)[elem_id] = 1;
      
      double low_fidelity_particle_distribution_function_value = 0.0;
      if (variance_reduction_parameters_.strategy == VarianceReductionParameters::Strategy::LocalMaxwellian)
      {
        mfem::Vector primitive_state(5);
        primitive_state(0) = number_density.at(particle.species)[elem_id];
        const mfem::Vector bulk_velocity_in_element(bulk_velocity.at(particle.species).GetColumn(elem_id), 3);
        primitive_state(0) = number_density.at(particle.species)(elem_id);
        primitive_state(1) = bulk_velocity_in_element(0);
        primitive_state(2) = bulk_velocity_in_element(1);
        primitive_state(3) = bulk_velocity_in_element(2);
        primitive_state(4) = temperature.at(particle.species)(elem_id);
        low_fidelity_particle_distribution_function_value = euler::evaluateMaxwellian(primitive_state,particle_velocity,particle.species);
      }
      else
        low_fidelity_particle_distribution_function_value = low_fidelity_operations.evaluateParticleDistributionFunction(low_fidelity_state,particle_position,particle_velocity,particle.element,low_fidelity_species_index);

      double noise_reducing_factor = (1 - low_fidelity_particle_distribution_function_value / particle.particle_distribution_function_value);
      velocity_in_element(0) += (particle.weight * particle.velocity(0) * noise_reducing_factor) / (variance_reduced_number_density * element_volume);
      velocity_in_element(1) += (particle.weight * particle.velocity(1) * noise_reducing_factor) / (variance_reduced_number_density * element_volume);
      velocity_in_element(2) += (particle.weight * particle.velocity(2) * noise_reducing_factor) / (variance_reduced_number_density * element_volume);
    }
    else
    {
      const double sum_weights = sum_of_weights_.at(particle.species)(elem_id);
      velocity_in_element.Add(particle.weight / sum_weights, particle.velocity);
    }
  }

  for (int elem_id = 0; elem_id < finite_element_space.GetNE(); ++elem_id)
  {
    const double element_volume = mesh.GetElementVolume(elem_id);
    for(int ispecies = 0; ispecies < low_fidelity_state.numSpecies(); ++ispecies)
    {
      const LowFidelitySpeciesState& current_species_state = low_fidelity_state.getSpeciesState(ispecies);
      Species current_species = current_species_state.getSpecies();
      mfem::Vector low_fidelity_integral_in_element(low_fidelity_integral.at(current_species).GetColumn(elem_id), 3);
      mfem::Vector velocity_in_element(variance_reduced_particle_bulk_velocity_.at(current_species).GetColumn(elem_id), 3);
      double variance_reduced_number_density = variance_reduced_particle_number_density_.at(current_species)(elem_id);
      if (variance_reduction_performed.at(current_species)[elem_id] == 1)
      {
        if (variance_reduction_parameters_.specified_lf_moments)
        {
          velocity_in_element(0) += variance_reduction_parameters_.reference_bulk_velocity[0];
          velocity_in_element(1) += variance_reduction_parameters_.reference_bulk_velocity[1];
          velocity_in_element(2) += variance_reduction_parameters_.reference_bulk_velocity[2];
        }
        else
        {
          velocity_in_element(0) += low_fidelity_integral_in_element(0) / (variance_reduced_number_density * element_volume);
          velocity_in_element(1) += low_fidelity_integral_in_element(1) / (variance_reduced_number_density * element_volume);
          velocity_in_element(2) += low_fidelity_integral_in_element(2) / (variance_reduced_number_density * element_volume);
        }
      }
    }
  } 

  return this->variance_reduced_particle_bulk_velocity_;
}

std::unordered_map<Species, mfem::Vector>& ParticleOperations::getTemperature(const ParticleContainer& particles, const bool sum_weights, const bool compute_bulk_velocity
) {

  for (auto & species_and_temperature : particle_temperature_)
    species_and_temperature.second = 0.0;

  if (sum_weights) this->sumParticleWeights_(particles);
  if (compute_bulk_velocity) this->getBulkVelocity(particles, false);

  for (const Particle& particle : particles) {
    if (not particle.is_alive) continue;

    const int elem_id = particle.element;
    const Species & species = particle.species;
    const double sum_weights = sum_of_weights_.at(species)(elem_id);
    if (sum_weights <= 0.0) continue;

    const mfem::Vector bulk_velocity_in_element(particle_bulk_velocity_.at(species).GetColumn(elem_id), 3);
    mfem::Vector fluctuation_velocity = particle.velocity;
    fluctuation_velocity -= bulk_velocity_in_element;

    const double norm_squared = fluctuation_velocity * fluctuation_velocity;

    // this only holds for constant weights
    const double number_of_samples = sum_weights / particle.weight;
    if (number_of_samples <= 1.) continue;

    const double bias_corrected_weight = particle.weight * number_of_samples / (number_of_samples - 1.);
    const double m_over_3kb = particle.species.mass / (velocity_dims_ * constants::boltzmann_constant);

    particle_temperature_.at(species)(elem_id) += norm_squared * bias_corrected_weight * m_over_3kb / sum_weights;
  }

  return this->particle_temperature_;
}

std::unordered_map<Species, double>& ParticleOperations::getAvgTemperature(const ParticleContainer& particles, const bool sum_weights, const bool compute_bulk_velocity
) {

  for (auto & species_and_temperature : avg_particle_temperature_)
    species_and_temperature.second = 0.0;

  if (sum_weights) this->sumParticleWeights_(particles);
  if (compute_bulk_velocity) this->getBulkVelocity(particles, false);

  for (const Particle& particle : particles) {
    if (not particle.is_alive) continue;

    const int elem_id = particle.element;
    const Species & species = particle.species;
    const double sum_weights = sum_of_weights_.at(species)(elem_id);
    if (sum_weights <= 0.0) continue;

    //const mfem::Vector bulk_velocity_in_element(particle_bulk_velocity_.at(species).GetColumn(elem_id), 3);
    mfem::Vector fluctuation_velocity = particle.velocity;
    //fluctuation_velocity -= bulk_velocity_in_element;

    const double norm_squared = fluctuation_velocity * fluctuation_velocity;

    // this only holds for constant weights
    const double number_of_samples = sum_weights / particle.weight;
    const double number_of_physical_particles = sum_of_weights_.at(particle.species).Sum();
    const double one_over_n_minus_one = 1./(number_of_physical_particles - 1);
    const double m_over_3kb = particle.species.mass / (velocity_dims_ * constants::boltzmann_constant);
    if (number_of_samples <= 1.) continue;

    avg_particle_temperature_.at(species) += one_over_n_minus_one * particle.weight * norm_squared * m_over_3kb;
  }

  return this->avg_particle_temperature_;
}

std::unordered_map<Species,mfem::Vector>& ParticleOperations::getVarianceReducedTemperature(
  const ParticleContainer& particles, 
  const LowFidelityState& low_fidelity_state,
  const DGEulerOperations& low_fidelity_operations
) {

  for (auto & species_and_temperature : variance_reduced_particle_temperature_)
    species_and_temperature.second = 0.0;

  mfem::FiniteElementSpace finite_element_space = discretization_.getFeSpace();
  mfem::Mesh &mesh = *finite_element_space.GetMesh();

  //TODO: Store and reuse from previous moments
  std::unordered_map<Species,mfem::Vector> low_fidelity_number_density_integral = low_fidelity_operations.integralForVarianceReducedNumberDensity(finite_element_space, low_fidelity_state);
  std::unordered_map<Species,mfem::DenseMatrix> low_fidelity_bulk_velocity_integral = low_fidelity_operations.integralForVarianceReducedBulkVelocity(finite_element_space, low_fidelity_state);
  std::unordered_map<Species,mfem::Vector> low_fidelity_temperature_integral = low_fidelity_operations.integralForVarianceReducedTemperature(finite_element_space, low_fidelity_state, velocity_dims_);

  this->sumParticleWeights_(particles);
  this->getVarianceReducedPostprocessors(particles,low_fidelity_state,low_fidelity_operations);

  std::unordered_map<Species, std::vector<char>> variance_reduction_performed; 
  std::unordered_map<Species, std::vector<char>> more_than_one_macro_particle;
  variance_reduction_performed.reserve(variance_reduced_particle_number_density_.size());
  more_than_one_macro_particle.reserve(variance_reduced_particle_number_density_.size());
  for (const auto& kv : variance_reduced_particle_number_density_)
  {
    variance_reduction_performed.emplace(kv.first, std::vector<char>(finite_element_space.GetNE(), 0));
    more_than_one_macro_particle.emplace(kv.first, std::vector<char>(finite_element_space.GetNE(), 1.0));
  }

  /////
  //Compute standard PIC moments
  std::unordered_map<Species, mfem::Vector> number_density     = this->getNumberDensity(particles);
  std::unordered_map<Species, mfem::DenseMatrix> bulk_velocity = this->getBulkVelocity(particles,true);
  std::unordered_map<Species, mfem::Vector> temperature        = this->getTemperature(particles,true,true);
  /////

  mfem::Vector particle_velocity(velocity_dims_);  
  for (const Particle& particle : particles) {
    if (not particle.is_alive) continue;

    for (int vel_dim = 0; vel_dim < velocity_dims_; ++vel_dim)
      particle_velocity(vel_dim) = particle.velocity(vel_dim);

    const int elem_id = particle.element;
    const double element_volume = mesh.GetElementVolume(elem_id);

    const mfem::Vector particle_position(particle.position.GetData(), dim_);

    double variance_reduced_number_density = variance_reduced_particle_number_density_.at(particle.species)(elem_id);
    mfem::Vector fluctuation_velocity = particle_velocity;
    for (int vel_dim = 0; vel_dim < velocity_dims_; ++vel_dim)
      fluctuation_velocity(vel_dim) = particle.velocity(vel_dim) - variance_reduced_particle_bulk_velocity_.at(particle.species)(vel_dim, elem_id);

    const double sum_weights = sum_of_weights_.at(particle.species)(elem_id);
    if (sum_weights <= 0.0) continue;
    const double number_of_macro_particles = sum_weights / particle.weight;

    const double norm_squared = fluctuation_velocity * fluctuation_velocity;
    const double m_over_3kb = particle.species.mass / (velocity_dims_ * constants::boltzmann_constant);
    const int low_fidelity_species_index =
      low_fidelity_state.getSpeciesIndex(particle.species);
    bool perform_variance_reduction = (low_fidelity_species_index >= 0);
    if (variance_reduction_parameters_.limit_variance_reduction)
    {
      perform_variance_reduction = (low_fidelity_species_index >= 0) &&
        (variance_reduced_postprocessors_.noise_reducing_factor
            .at(particle.species)(elem_id) < 1.0);
    }

    if (number_of_macro_particles > 1)
    {
      if (perform_variance_reduction)
      {
        variance_reduction_performed.at(particle.species)[elem_id] = 1;

        double low_fidelity_particle_distribution_function_value = 0.0;
        if (variance_reduction_parameters_.strategy == VarianceReductionParameters::Strategy::LocalMaxwellian)
        {
          mfem::Vector primitive_state(5);
          const mfem::Vector bulk_velocity_in_element(bulk_velocity.at(particle.species).GetColumn(elem_id), 3);
          primitive_state(0) = number_density.at(particle.species)[elem_id];
          primitive_state(1) = bulk_velocity_in_element(0);
          primitive_state(2) = bulk_velocity_in_element(1);
          primitive_state(3) = bulk_velocity_in_element(2);
          primitive_state(4) = temperature.at(particle.species)(elem_id);
          low_fidelity_particle_distribution_function_value = euler::evaluateMaxwellian(primitive_state,particle_velocity,particle.species);
        }
        else
          low_fidelity_particle_distribution_function_value = low_fidelity_operations.evaluateParticleDistributionFunction(low_fidelity_state,particle_position,particle_velocity,particle.element,low_fidelity_species_index);

        double noise_reducing_factor = (1 - low_fidelity_particle_distribution_function_value / particle.particle_distribution_function_value);
        variance_reduced_particle_temperature_.at(particle.species)(elem_id) += number_of_macro_particles / (number_of_macro_particles - 1) * m_over_3kb * norm_squared * particle.weight * noise_reducing_factor / (variance_reduced_number_density * element_volume) ;
      }
      else
      {
        const double bias_corrected_weight = particle.weight * number_of_macro_particles / (number_of_macro_particles - 1.);
        variance_reduced_particle_temperature_.at(particle.species)(elem_id) += norm_squared * bias_corrected_weight * m_over_3kb / sum_weights;
      }
    }
    else
    {
      more_than_one_macro_particle.at(particle.species)[elem_id] = 0.0;
    }
  }

  for (int elem_id = 0; elem_id < finite_element_space.GetNE(); ++elem_id)
  {
    const double element_volume = mesh.GetElementVolume(elem_id);
    for(int ispecies = 0; ispecies < low_fidelity_state.numSpecies(); ++ispecies)
    {
      const LowFidelitySpeciesState& current_species_state = low_fidelity_state.getSpeciesState(ispecies);
      Species current_species = current_species_state.getSpecies();
      if (variance_reduction_performed.at(current_species)[elem_id] == 1) 
      {
        if (variance_reduction_parameters_.specified_lf_moments)
        {
          variance_reduced_particle_temperature_.at(current_species)(elem_id) += variance_reduction_parameters_.reference_temperature;
        }
        else
        {
          const double m_over_3kb = current_species.mass / (velocity_dims_ * constants::boltzmann_constant);
          double number_density = variance_reduced_particle_number_density_.at(current_species)(elem_id);

          double x_bulk_velocity = variance_reduced_particle_bulk_velocity_.at(current_species)(0,elem_id);
          double y_bulk_velocity = variance_reduced_particle_bulk_velocity_.at(current_species)(1,elem_id);
          double z_bulk_velocity = variance_reduced_particle_bulk_velocity_.at(current_species)(2,elem_id);
          double bulk_velocity_mag_squared = x_bulk_velocity * x_bulk_velocity + y_bulk_velocity * y_bulk_velocity + z_bulk_velocity * z_bulk_velocity;
          const double bulk_velocity_dot_low_fidelity_bulk_velocity_integral 
          = x_bulk_velocity * low_fidelity_bulk_velocity_integral.at(current_species)(0,elem_id)
          + y_bulk_velocity * low_fidelity_bulk_velocity_integral.at(current_species)(1,elem_id)
          + z_bulk_velocity * low_fidelity_bulk_velocity_integral.at(current_species)(2,elem_id);

          variance_reduced_particle_temperature_.at(current_species)(elem_id) 
            += more_than_one_macro_particle.at(current_species)[elem_id] * 
            (low_fidelity_temperature_integral.at(current_species)(elem_id) / (number_density * element_volume)
            + m_over_3kb * bulk_velocity_mag_squared * low_fidelity_number_density_integral.at(current_species)(elem_id) / (number_density * element_volume)
            - m_over_3kb * 2.0 * bulk_velocity_dot_low_fidelity_bulk_velocity_integral / (number_density * element_volume));
        }
      }
    }
  } 
  return this->variance_reduced_particle_temperature_;
}

void ParticleOperations::sumParticleWeights_(
  const ParticleContainer& particles
) {
    for (auto & species_and_sum_of_weights : sum_of_weights_)
      species_and_sum_of_weights.second = 0.0;
    for (const Particle& particle : particles) {
      if (not particle.is_alive) continue;

      const int elem_id = particle.element;
      const Species & species = particle.species;
      sum_of_weights_.at(species)(elem_id) += particle.weight;
    }
}


VarianceReducedPostprocessors ParticleOperations::getVarianceReducedPostprocessors(
  const ParticleContainer& particles, 
  const LowFidelityState& low_fidelity_state,
  const DGEulerOperations& low_fidelity_operations
) {

  for (auto & species_and_number_density : variance_reduced_postprocessors_.noise_reducing_factor)
    species_and_number_density.second = 0.0;

  for (auto & species_and_number_density : variance_reduced_postprocessors_.particle_distribution_function)
    species_and_number_density.second = std::numeric_limits<double>::infinity();

  for (auto & species_and_number_density : variance_reduced_postprocessors_.low_fidelity_particle_distribution_function)
    species_and_number_density.second = std::numeric_limits<double>::infinity();

  mfem::FiniteElementSpace finite_element_space = discretization_.getFeSpace();

  if (variance_reduction_parameters_.strategy == VarianceReductionParameters::Strategy::LocalMaxwellian)
  {
    this->getNumberDensity(particles);
    this->getBulkVelocity(particles,true);
    this->getTemperature(particles,true,true);
  }

  mfem::Vector particle_velocity(velocity_dims_);  
  for (const Particle& particle : particles) {
    if (not particle.is_alive) continue;

    for (int vel_dim = 0; vel_dim < velocity_dims_; ++vel_dim)
      particle_velocity(vel_dim) = particle.velocity(vel_dim);

    const int elem_id = particle.element;
    const mfem::Vector particle_position(particle.position.GetData(), dim_);
    const int low_fidelity_species_index =
      low_fidelity_state.getSpeciesIndex(particle.species);
    double noise_reducing_factor = -10;
    double low_fidelity_particle_distribution_function_value = 0.0;
    if (low_fidelity_species_index >= 0)
    {
      if (variance_reduction_parameters_.strategy == VarianceReductionParameters::Strategy::LocalMaxwellian)
      {
        mfem::Vector primitive_state(5);
        primitive_state(0) = particle_number_density_.at(particle.species)[elem_id];
        const mfem::Vector bulk_velocity_in_element(particle_bulk_velocity_.at(particle.species).GetColumn(elem_id), 3);
        primitive_state(0) = particle_number_density_.at(particle.species)(elem_id);
        primitive_state(1) = bulk_velocity_in_element(0);
        primitive_state(2) = bulk_velocity_in_element(1);
        primitive_state(3) = bulk_velocity_in_element(2);
        primitive_state(4) = particle_temperature_.at(particle.species)(elem_id);
        low_fidelity_particle_distribution_function_value = euler::evaluateMaxwellian(primitive_state,particle_velocity,particle.species);
      }
      else
        low_fidelity_particle_distribution_function_value = low_fidelity_operations.evaluateParticleDistributionFunction(low_fidelity_state,particle_position,particle_velocity,particle.element,low_fidelity_species_index);
      noise_reducing_factor = (1 - low_fidelity_particle_distribution_function_value / particle.particle_distribution_function_value);
    }
    variance_reduced_postprocessors_.noise_reducing_factor.at(particle.species)(elem_id) = std::max(
      std::abs(noise_reducing_factor),
      variance_reduced_postprocessors_.noise_reducing_factor.at(particle.species)(elem_id));

    variance_reduced_postprocessors_.particle_distribution_function.at(particle.species)(elem_id) = std::min(
      particle.particle_distribution_function_value,
      variance_reduced_postprocessors_.particle_distribution_function.at(particle.species)(elem_id));

    variance_reduced_postprocessors_.low_fidelity_particle_distribution_function.at(particle.species)(elem_id) = std::min(
      low_fidelity_particle_distribution_function_value,
      variance_reduced_postprocessors_.low_fidelity_particle_distribution_function.at(particle.species)(elem_id));
  }
  return this->variance_reduced_postprocessors_;
}

void ParticleOperations::updateParticleDistributionFunctionValue(
  ParticleContainer& particles,
  const ParticleContainer& reference_particles)
{
  std::unordered_map<Species, mfem::Vector> number_density     = this->getNumberDensity(reference_particles);
  std::unordered_map<Species, mfem::DenseMatrix> bulk_velocity = this->getBulkVelocity(reference_particles,true);
  std::unordered_map<Species, mfem::Vector> temperature        = this->getTemperature(reference_particles,false,false);

  mfem::Vector particle_velocity(velocity_dims_);  
  for (Particle& particle : particles) {
    if (not particle.is_alive) continue;

    for (int vel_dim = 0; vel_dim < velocity_dims_; ++vel_dim)
      particle_velocity(vel_dim) = particle.velocity(vel_dim);

    const int elem_id = particle.element;
    const Species & species = particle.species;
    const mfem::Vector bulk_velocity_in_element(bulk_velocity.at(species).GetColumn(elem_id), 3);

    mfem::Vector primitive_state(5);
    primitive_state(0) = number_density.at(species)(elem_id);
    primitive_state(1) = bulk_velocity_in_element(0);
    primitive_state(2) = bulk_velocity_in_element(1);
    primitive_state(3) = bulk_velocity_in_element(2);
    primitive_state(4) = temperature.at(species)(elem_id);
    double particle_distribution_function_value = euler::evaluateMaxwellian(primitive_state,particle_velocity,species);
    particle.particle_distribution_function_value = particle_distribution_function_value;
  }
}

void ParticleOperations::updateParticleDistributionFunctionValue(
  ParticleContainer& particles,
  double number_density,
  mfem::Vector bulk_velocity,
  double temperature)
{
  mfem::Vector particle_velocity(velocity_dims_);  
  for (Particle& particle : particles) {
    if (not particle.is_alive) continue;

    for (int vel_dim = 0; vel_dim < velocity_dims_; ++vel_dim)
      particle_velocity(vel_dim) = particle.velocity(vel_dim);

    const Species & species = particle.species;
    mfem::Vector primitive_state(5);
    primitive_state(0) = number_density;
    primitive_state(1) = bulk_velocity(0);
    primitive_state(2) = bulk_velocity(1);
    primitive_state(3) = bulk_velocity(2);
    primitive_state(4) = temperature;
    double particle_distribution_function_value = euler::evaluateMaxwellian(primitive_state,particle_velocity,species);
    particle.particle_distribution_function_value = particle_distribution_function_value;
  }
}

std::unordered_map<Species, double>& ParticleOperations::getAvgVarianceReducedNumberDensity(
  const ParticleContainer& particles, 
  const LowFidelityState& low_fidelity_state,
  const DGEulerOperations& low_fidelity_operations
) {

  for (auto & species_and_number_density : avg_variance_reduced_particle_number_density_)
    species_and_number_density.second = 0.0;

  mfem::FiniteElementSpace finite_element_space = discretization_.getFeSpace();
  mfem::Mesh &mesh = *finite_element_space.GetMesh();

  std::unordered_map<Species,mfem::Vector> low_fidelity_integral = low_fidelity_operations.integralForVarianceReducedNumberDensity(finite_element_space, low_fidelity_state);
  this->getVarianceReducedPostprocessors(particles,low_fidelity_state,low_fidelity_operations);

  for (const Particle& particle : particles) {
    if (not particle.is_alive) continue;

    const mfem::Vector particle_position(particle.position.GetData(), dim_);
    const int low_fidelity_species_index =
      low_fidelity_state.getSpeciesIndex(particle.species);
    double low_fidelity_particle_distribution_function_value = low_fidelity_operations.evaluateParticleDistributionFunction(low_fidelity_state,particle_position,particle.velocity,particle.element,low_fidelity_species_index);
    double noise_reducing_factor = (1 - low_fidelity_particle_distribution_function_value / particle.particle_distribution_function_value);
    variance_reduced_particle_number_density_.at(particle.species) += (particle.weight * noise_reducing_factor) / getMeshVolume(mesh);
  }

  for (int elem_id = 0; elem_id < finite_element_space.GetNE(); ++elem_id) {
    const double element_volume = mesh.GetElementVolume(elem_id);
    for (int ispecies = 0; ispecies < low_fidelity_state.numSpecies(); ++ispecies) {
    const Species species = low_fidelity_state.getSpeciesState(ispecies).getSpecies();
      avg_variance_reduced_particle_number_density_.at(species) +=
        low_fidelity_integral.at(species)(elem_id) / element_volume;
    }
  }
  return this->avg_variance_reduced_particle_number_density_;
}

std::unordered_map<Species, mfem::Vector>& ParticleOperations::getAvgVarianceReducedBulkVelocity(
  const ParticleContainer& particles, 
  const LowFidelityState& low_fidelity_state,
  const DGEulerOperations& low_fidelity_operations
) {

  for (auto & species_and_bulk_velocity : avg_variance_reduced_particle_bulk_velocity_)
    species_and_bulk_velocity.second = 0.0;
  mfem::FiniteElementSpace finite_element_space = discretization_.getFeSpace();
  std::unordered_map<Species, mfem::DenseMatrix> low_fidelity_integral = low_fidelity_operations.integralForVarianceReducedBulkVelocity(finite_element_space, low_fidelity_state);
  this->sumParticleWeights_(particles);

  for (const Particle& particle : particles) {
    if (not particle.is_alive) continue;

    const double number_of_physical_particles = sum_of_weights_.at(particle.species).Sum();

    const mfem::Vector particle_position(particle.position.GetData(), dim_);
    const int low_fidelity_species_index =
      low_fidelity_state.getSpeciesIndex(particle.species);
    mfem::Vector velocity_in_element(avg_variance_reduced_particle_bulk_velocity_.at(particle.species));
    double low_fidelity_particle_distribution_function_value = low_fidelity_operations.evaluateParticleDistributionFunction(low_fidelity_state,particle_position,particle.velocity,particle.element,low_fidelity_species_index);
    double noise_reducing_factor = (1 - low_fidelity_particle_distribution_function_value / particle.particle_distribution_function_value);
    avg_variance_reduced_particle_bulk_velocity_.at(particle.species)(0) += (particle.weight * particle.velocity(0) * noise_reducing_factor) / number_of_physical_particles;
    avg_variance_reduced_particle_bulk_velocity_.at(particle.species)(1) += (particle.weight * particle.velocity(1) * noise_reducing_factor) / number_of_physical_particles; 
    avg_variance_reduced_particle_bulk_velocity_.at(particle.species)(2) += (particle.weight * particle.velocity(2) * noise_reducing_factor) / number_of_physical_particles;
  }

  for (int elem_id = 0; elem_id < finite_element_space.GetNE(); ++elem_id)
  {
    for(int ispecies = 0; ispecies < low_fidelity_state.numSpecies(); ++ispecies)
    {
      const LowFidelitySpeciesState& current_species_state = low_fidelity_state.getSpeciesState(ispecies);
      Species current_species = current_species_state.getSpecies();
      const double number_of_physical_particles = sum_of_weights_.at(current_species)(0);
      mfem::Vector low_fidelity_integral_in_element(low_fidelity_integral.at(current_species).GetColumn(elem_id), 3);
      mfem::Vector velocity_in_element(avg_variance_reduced_particle_bulk_velocity_.at(current_species));
      avg_variance_reduced_particle_bulk_velocity_.at(current_species)(0) += low_fidelity_integral_in_element(0) / number_of_physical_particles;
      avg_variance_reduced_particle_bulk_velocity_.at(current_species)(1) += low_fidelity_integral_in_element(1) / number_of_physical_particles;
      avg_variance_reduced_particle_bulk_velocity_.at(current_species)(2) += low_fidelity_integral_in_element(2) / number_of_physical_particles;
    }
  } 

  return this->avg_variance_reduced_particle_bulk_velocity_;
}

std::unordered_map<Species,double>& ParticleOperations::getAvgVarianceReducedTemperature(
  const ParticleContainer& particles, 
  const LowFidelityState& low_fidelity_state,
  const DGEulerOperations& low_fidelity_operations
) {

  for (auto & species_and_temperature : avg_variance_reduced_particle_temperature_)
    species_and_temperature.second = 0.0;

  mfem::FiniteElementSpace finite_element_space = discretization_.getFeSpace();
  //mfem::Mesh &mesh = *finite_element_space.GetMesh();

  std::unordered_map<Species,mfem::Vector> low_fidelity_number_density_integral = low_fidelity_operations.integralForVarianceReducedNumberDensity(finite_element_space, low_fidelity_state);
  std::unordered_map<Species,mfem::DenseMatrix> low_fidelity_bulk_velocity_integral = low_fidelity_operations.integralForVarianceReducedBulkVelocity(finite_element_space, low_fidelity_state);
  std::unordered_map<Species,mfem::Vector> low_fidelity_temperature_integral = low_fidelity_operations.integralForVarianceReducedTemperature(finite_element_space, low_fidelity_state, velocity_dims_);

  this->sumParticleWeights_(particles);
  for (const Particle& particle : particles) {
    if (not particle.is_alive) continue;

    //const int elem_id = particle.element;
    const mfem::Vector particle_position(particle.position.GetData(), dim_);
    //const double sum_weights = sum_of_weights_.at(particle.species)(elem_id);
    //const double number_of_macro_particles = sum_weights / particle.weight;
    const double number_of_physical_particles = sum_of_weights_.at(particle.species).Sum();
    const double one_over_n_minus_one = 1./(number_of_physical_particles - 1);

    //could replace with zero here if known
    // double x_bulk_velocity = particle_bulk_velocity_.at(particle.species)(0,elem_id);
    // double y_bulk_velocity = particle_bulk_velocity_.at(particle.species)(1,elem_id);
    // double z_bulk_velocity = particle_bulk_velocity_.at(particle.species)(2,elem_id);
    double x_bulk_velocity = 0.0;
    double y_bulk_velocity = 0.0;
    double z_bulk_velocity = 0.0;

    mfem::Vector fluctuation_velocity = particle.velocity;
    fluctuation_velocity(0) -= x_bulk_velocity;
    fluctuation_velocity(1) -= y_bulk_velocity;
    fluctuation_velocity(2) -= z_bulk_velocity;
    const double norm_squared = fluctuation_velocity * fluctuation_velocity;
    const double m_over_3kb = particle.species.mass / (3.0 * constants::boltzmann_constant);
    const int low_fidelity_species_index =
      low_fidelity_state.getSpeciesIndex(particle.species);
    double low_fidelity_particle_distribution_function_value = low_fidelity_operations.evaluateParticleDistributionFunction(low_fidelity_state,particle_position,particle.velocity,particle.element,low_fidelity_species_index);
    double noise_reducing_factor = (1 - low_fidelity_particle_distribution_function_value / particle.particle_distribution_function_value);
    //variance_reduced_particle_temperature_.at(particle.species)(elem_id) += number_of_macro_particles / (number_of_macro_particles - 1) * m_over_3kb * norm_squared * particle.weight * noise_reducing_factor / number_of_physical_particles;
    avg_variance_reduced_particle_temperature_.at(particle.species) += one_over_n_minus_one * m_over_3kb * norm_squared * particle.weight * noise_reducing_factor;
  }

  for (int elem_id = 0; elem_id < finite_element_space.GetNE(); ++elem_id)
  {
    for(int ispecies = 0; ispecies < low_fidelity_state.numSpecies(); ++ispecies)
    {
      const LowFidelitySpeciesState& current_species_state = low_fidelity_state.getSpeciesState(ispecies);
      Species current_species = current_species_state.getSpecies();
      const double number_of_physical_particles = sum_of_weights_.at(current_species).Sum();
      const double m_over_3kb = current_species.mass / (3.0 * constants::boltzmann_constant);

      double x_bulk_velocity = avg_variance_reduced_particle_bulk_velocity_.at(current_species)(0);
      double y_bulk_velocity = avg_variance_reduced_particle_bulk_velocity_.at(current_species)(1);
      double z_bulk_velocity = avg_variance_reduced_particle_bulk_velocity_.at(current_species)(2);
      double bulk_velocity_mag_squared = x_bulk_velocity * x_bulk_velocity + y_bulk_velocity * y_bulk_velocity + z_bulk_velocity * z_bulk_velocity;
      const double bulk_velocity_dot_low_fidelity_bulk_velocity_integral 
      = x_bulk_velocity * low_fidelity_bulk_velocity_integral.at(current_species)(0,elem_id)
      + y_bulk_velocity * low_fidelity_bulk_velocity_integral.at(current_species)(1,elem_id)
      + z_bulk_velocity * low_fidelity_bulk_velocity_integral.at(current_species)(2,elem_id);

      avg_variance_reduced_particle_temperature_.at(current_species)
        += 
        (low_fidelity_temperature_integral.at(current_species)(elem_id) / number_of_physical_particles
        + m_over_3kb * bulk_velocity_mag_squared * low_fidelity_number_density_integral.at(current_species)(elem_id) / number_of_physical_particles
        - m_over_3kb * 2.0 * bulk_velocity_dot_low_fidelity_bulk_velocity_integral / number_of_physical_particles);
    }
  } 
  return this->avg_variance_reduced_particle_temperature_;
}

void ParticleOperations::updateParticleDistributionFunctionValue(
    ParticleContainer& particles,
    const VelocityHistogram& hist)
{
    for (Particle& particle : particles) {

        if (!particle.is_alive) continue;

        particle.particle_distribution_function_value =
            hist.evaluate(particle.velocity(0));
    }
}

} // namespace mfpic
