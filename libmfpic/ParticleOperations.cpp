#include <libmfpic/BuildVarianceReductionParametersFromYaml.hpp>
#include <libmfpic/DGEulerOperations.hpp>
#include <libmfpic/Constants.hpp>
#include <libmfpic/Errors.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/IntegratedCharge.hpp>
#include <libmfpic/LowFidelityOperations.hpp>
#include <libmfpic/LowFidelityState.hpp>
#include <libmfpic/MeshUtilities.hpp>
#include <libmfpic/ParticleContainer.hpp>
#include <libmfpic/ParticleOperations.hpp>
#include <libmfpic/PeriodicParticleBoundary.hpp>
#include <libmfpic/Species.hpp>

#include <mfem/mfem.hpp>

#include <algorithm>
#include <limits>
#include <random>
#include <ranges>
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
    particle_moments_.number_density.insert({species, mfem::Vector(mesh.GetNE())});
    particle_moments_.bulk_velocity.insert({species, mfem::DenseMatrix(3, mesh.GetNE())});
    particle_moments_.temperature.insert({species, mfem::Vector(mesh.GetNE())});
    variance_reduced_particle_moments_.number_density.insert({species, mfem::Vector(mesh.GetNE())});
    variance_reduced_particle_moments_.bulk_velocity.insert({species, mfem::DenseMatrix(3, mesh.GetNE())});
    variance_reduced_particle_moments_.temperature.insert({species, mfem::Vector(mesh.GetNE())});
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
  return assembleVarianceReducedCharge_(current_particles, &low_fidelity_state, &low_fidelity_operations);
}

IntegratedCharge ParticleOperations::assembleVarianceReducedCharge(
  const ParticleContainer& current_particles
) const {
  return assembleVarianceReducedCharge_(current_particles, nullptr, nullptr);
}

IntegratedCharge ParticleOperations::assembleVarianceReducedCharge_(
  const ParticleContainer& current_particles,
  const LowFidelityState* low_fidelity_state,
  const LowFidelityOperations* low_fidelity_operations
) const {
  mfem::IntegrationPoint integration_point;
  mfem::Array<int> vector_dofs;

  mfem::FiniteElementSpace& finite_element_space = discretization_.getFeSpace();
  mfem::Mesh& mesh = *finite_element_space.GetMesh();

  const std::unordered_map<Species, std::vector<char>> variance_reduction_elements =
    computeVarianceReductionElements_(current_particles, low_fidelity_state, low_fidelity_operations);

  // A DOF uses variance reduction only if every element it touches does, so that the particle and
  // control variate contributions to it cover the same region.
  std::unordered_map<Species, std::vector<char>> variance_reduction_dofs;
  for (const auto& [species, elements] : variance_reduction_elements) {
    std::vector<char> dofs(finite_element_space.GetNDofs(), 1);
    for (int element_id = 0; element_id < finite_element_space.GetNE(); ++element_id) {
      if (elements[element_id]) continue;
      finite_element_space.GetElementVDofs(element_id, vector_dofs);
      for (const int global_dof : vector_dofs) dofs[global_dof] = 0;
    }
    variance_reduction_dofs.emplace(species, std::move(dofs));
  }

  IntegratedCharge integrated_charge(discretization_);
  integrated_charge.setIntegratedChargeValue(0.0);

  for (const Particle& particle : current_particles) {
    if (!particle.is_alive) continue;

    const int element_id = particle.element;
    const Species& particle_species = particle.species;
    const double particle_charge = particle_species.charge;
    const std::vector<char>& dofs_using_variance_reduction = variance_reduction_dofs.at(particle_species);

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

    const bool any_dof_uses_variance_reduction = std::ranges::any_of(
      vector_dofs, [&](const int global_dof) { return dofs_using_variance_reduction[global_dof] == 1; });
    const double noise_reducing_factor = any_dof_uses_variance_reduction
      ? 1.0 - evaluateControlVariate_(particle, low_fidelity_state, low_fidelity_operations) / particle.particle_distribution_function_value
      : 1.0;

    for (int local_dof = 0; local_dof < finite_element->GetDof(); ++local_dof) {
      const int global_dof = vector_dofs[local_dof];
      const double factor = dofs_using_variance_reduction[global_dof] == 1 ? noise_reducing_factor : 1.0;
      integrated_charge.addIntegratedChargeValue(global_dof, particle.weight * particle_charge * shape_functions(local_dof) * factor);
    }
  }

  for (const auto& [species, dofs] : variance_reduction_dofs) {
    if (!hasControlVariate_(species, low_fidelity_state)) continue;
    const mfem::Vector control_variate_charge =
      assembleControlVariateCharge_(species, low_fidelity_state, low_fidelity_operations);
    for (int global_dof = 0; global_dof < finite_element_space.GetNDofs(); ++global_dof) {
      if (dofs[global_dof] == 1) {
        integrated_charge.addIntegratedChargeValue(global_dof, control_variate_charge(global_dof));
      }
    }
  }
  return integrated_charge;
}

ParticleMoments& ParticleOperations::getParticleMoments(const ParticleContainer& particles
) {
  sumParticleWeights_(particles);
  getNumberDensity(particles);
  getBulkVelocity(particles,false);
  getTemperature(particles,false,false);
  return particle_moments_;
}

ParticleMoments& ParticleOperations::getVarianceReducedParticleMoments(
  const ParticleContainer& particles,
  const LowFidelityState& low_fidelity_state,
  const DGEulerOperations& low_fidelity_operations,
  const bool compute_standard_pic_moments
) {
  return getVarianceReducedParticleMoments_(particles, &low_fidelity_state, &low_fidelity_operations, compute_standard_pic_moments);
}

ParticleMoments& ParticleOperations::getVarianceReducedParticleMoments(
  const ParticleContainer& particles,
  const bool compute_standard_pic_moments
) {
  return getVarianceReducedParticleMoments_(particles, nullptr, nullptr, compute_standard_pic_moments);
}

ParticleMoments& ParticleOperations::getVarianceReducedParticleMoments_(
  const ParticleContainer& particles,
  const LowFidelityState* low_fidelity_state,
  const DGEulerOperations* low_fidelity_operations,
  const bool compute_standard_pic_moments
) {
  if (compute_standard_pic_moments)
    getParticleMoments(particles);

  variance_reduction_elements_ = computeVarianceReductionElements_(particles, low_fidelity_state, low_fidelity_operations);
  control_variate_moment_integrals_ = computeControlVariateMomentIntegrals_(low_fidelity_state, low_fidelity_operations);

  getVarianceReducedNumberDensity(particles,low_fidelity_state,low_fidelity_operations);
  getVarianceReducedBulkVelocity(particles,low_fidelity_state,low_fidelity_operations);
  getVarianceReducedTemperature(particles,low_fidelity_state,low_fidelity_operations);
  return variance_reduced_particle_moments_;
}

std::unordered_map<Species, mfem::Vector>& ParticleOperations::getNumberDensity(const ParticleContainer& particles
) {

  for (auto & species_and_number_density : particle_moments_.number_density)
    species_and_number_density.second = 0.0;

  mfem::Array<int> vector_dofs;
  mfem::FiniteElementSpace finite_element_space = discretization_.getFeSpace();
  mfem::Mesh &mesh = *finite_element_space.GetMesh();

  for (const Particle& particle : particles) {
    if (not particle.is_alive) continue;

    const int elem_id = particle.element;
    const Species & species = particle.species;
    particle_moments_.number_density.at(species)(elem_id) += particle.weight / mesh.GetElementVolume(elem_id);
  }

  return particle_moments_.number_density;
}

std::unordered_map<Species, mfem::Vector>& ParticleOperations::getVarianceReducedNumberDensity(
  const ParticleContainer& particles,
  const LowFidelityState* low_fidelity_state,
  const LowFidelityOperations* low_fidelity_operations
) {

  for (auto & species_and_number_density : variance_reduced_particle_moments_.number_density)
    species_and_number_density.second = 0.0;

  mfem::FiniteElementSpace& finite_element_space = discretization_.getFeSpace();
  mfem::Mesh &mesh = *finite_element_space.GetMesh();

  for (const Particle& particle : particles) {
    if (not particle.is_alive) continue;

    const int elem_id = particle.element;
    const double element_volume = mesh.GetElementVolume(elem_id);

    if (variance_reduction_elements_.at(particle.species)[elem_id])
    {
      double low_fidelity_particle_distribution_function_value = evaluateControlVariate_(particle,low_fidelity_state,low_fidelity_operations);
      double noise_reducing_factor = (1 - low_fidelity_particle_distribution_function_value / particle.particle_distribution_function_value);
      variance_reduced_particle_moments_.number_density.at(particle.species)(elem_id) += (particle.weight * noise_reducing_factor) / element_volume;
    }
    else
    {
      variance_reduced_particle_moments_.number_density.at(particle.species)(elem_id) += particle.weight / element_volume;
    }
  }
  for (const auto& [species, elements] : variance_reduction_elements_) {
    for (int elem_id = 0; elem_id < finite_element_space.GetNE(); ++elem_id) {
      if (elements[elem_id]) {
        variance_reduced_particle_moments_.number_density.at(species)(elem_id) +=
          control_variate_moment_integrals_.number_density.at(species)(elem_id) / mesh.GetElementVolume(elem_id);
      }
    }
  }
  return this->variance_reduced_particle_moments_.number_density;
}

std::unordered_map<Species, mfem::DenseMatrix>& ParticleOperations::getBulkVelocity(const ParticleContainer& particles, const bool sum_weights
) {

  for (auto & species_and_bulk_velocity : particle_moments_.bulk_velocity)
    species_and_bulk_velocity.second = 0.0;

  if (sum_weights) this->sumParticleWeights_(particles);

  for (const Particle& particle : particles) {
    if (not particle.is_alive) continue;

    const int elem_id = particle.element;
    const Species & species = particle.species;
    const double sum_weights = sum_of_weights_.at(species)(elem_id);
    if (sum_weights <= 0.0) continue;

    mfem::Vector velocity_in_element(particle_moments_.bulk_velocity.at(species).GetColumn(elem_id), 3);
    velocity_in_element.Add(particle.weight / sum_weights, particle.velocity);
  }

  return this->particle_moments_.bulk_velocity;
}

std::unordered_map<Species, mfem::DenseMatrix>& ParticleOperations::getVarianceReducedBulkVelocity(
  const ParticleContainer& particles,
  const LowFidelityState* low_fidelity_state,
  const LowFidelityOperations* low_fidelity_operations
) {

  for (auto & species_and_bulk_velocity : variance_reduced_particle_moments_.bulk_velocity)
    species_and_bulk_velocity.second = 0.0;

  this->sumParticleWeights_(particles);
  mfem::FiniteElementSpace& finite_element_space = discretization_.getFeSpace();
  mfem::Mesh &mesh = *finite_element_space.GetMesh();

  auto uses_variance_reduced_velocity = [&](const Species& species, const int elem_id) {
    const double number_density = variance_reduced_particle_moments_.number_density.at(species)(elem_id);
    return variance_reduction_elements_.at(species)[elem_id] && number_density > 0.0 && std::isfinite(number_density);
  };

  for (const Particle& particle : particles) {
    if (not particle.is_alive) continue;

    const int elem_id = particle.element;
    const double element_volume = mesh.GetElementVolume(elem_id);

    mfem::Vector velocity_in_element(variance_reduced_particle_moments_.bulk_velocity.at(particle.species).GetColumn(elem_id), 3);
    double variance_reduced_number_density = variance_reduced_particle_moments_.number_density.at(particle.species)(elem_id);
    if (uses_variance_reduced_velocity(particle.species, elem_id))
    {
      double low_fidelity_particle_distribution_function_value = evaluateControlVariate_(particle,low_fidelity_state,low_fidelity_operations);
      double noise_reducing_factor = (1 - low_fidelity_particle_distribution_function_value / particle.particle_distribution_function_value);
      for (int vel_dim = 0; vel_dim < velocity_dims_; ++vel_dim)
        velocity_in_element(vel_dim) += particle.weight * particle.velocity(vel_dim) * noise_reducing_factor / (variance_reduced_number_density * element_volume);
    }
    else
    {
      const double sum_weights = sum_of_weights_.at(particle.species)(elem_id);
      velocity_in_element.Add(particle.weight / sum_weights, particle.velocity);
    }
  }

  for (const Species& species : std::views::keys(variance_reduction_elements_))
  {
    for (int elem_id = 0; elem_id < finite_element_space.GetNE(); ++elem_id)
    {
      if (!uses_variance_reduced_velocity(species, elem_id)) continue;

      const double element_volume = mesh.GetElementVolume(elem_id);
      const mfem::DenseMatrix& momentum_integral = control_variate_moment_integrals_.momentum.at(species);
      mfem::Vector velocity_in_element(variance_reduced_particle_moments_.bulk_velocity.at(species).GetColumn(elem_id), 3);
      double variance_reduced_number_density = variance_reduced_particle_moments_.number_density.at(species)(elem_id);
      for (int vel_dim = 0; vel_dim < velocity_dims_; ++vel_dim)
        velocity_in_element(vel_dim) += momentum_integral(vel_dim, elem_id) / (variance_reduced_number_density * element_volume);
    }
  }

  return this->variance_reduced_particle_moments_.bulk_velocity;
}

std::unordered_map<Species, mfem::Vector>& ParticleOperations::getTemperature(const ParticleContainer& particles, const bool sum_weights, const bool compute_bulk_velocity
) {

  for (auto & species_and_temperature : particle_moments_.temperature)
    species_and_temperature.second = 0.0;

  if (sum_weights) this->sumParticleWeights_(particles);
  if (compute_bulk_velocity) this->getBulkVelocity(particles, false);

  std::unordered_map<Species, mfem::Vector> sum_of_squared_weights = sum_of_weights_;
  for (mfem::Vector& sum_of_squared_weights_for_species : std::views::values(sum_of_squared_weights)) {
    sum_of_squared_weights_for_species = 0.0;
  }
  for (const Particle& particle : particles) {
    if (not particle.is_alive) continue;

    sum_of_squared_weights.at(particle.species)(particle.element) += particle.weight * particle.weight;
  }

  for (const Particle& particle : particles) {
    if (not particle.is_alive) continue;

    const int elem_id = particle.element;
    const Species & species = particle.species;
    const double sum_of_weights_in_element = sum_of_weights_.at(species)(elem_id);
    if (sum_of_weights_in_element <= 0.0) continue;

    const mfem::Vector bulk_velocity_in_element(particle_moments_.bulk_velocity.at(species).GetColumn(elem_id), velocity_dims_);
    mfem::Vector fluctuation_velocity(velocity_dims_);
    subtract(mfem::Vector(particle.velocity.GetData(), velocity_dims_), bulk_velocity_in_element, fluctuation_velocity);
    const double norm_squared = fluctuation_velocity * fluctuation_velocity;

    const double sum_of_squared_weights_in_element = sum_of_squared_weights.at(species)(elem_id);
    const double sum_of_weights_in_element_squared = std::pow(sum_of_weights_in_element, 2.0);
    const double effective_num_particles = sum_of_weights_in_element_squared / sum_of_squared_weights_in_element;
    if (effective_num_particles == 1.0) {
      particle_moments_.temperature.at(species)(elem_id) = 0.0;
    }
    else {
      const double bias_corrected_weight = effective_num_particles / (effective_num_particles - 1.0) * particle.weight;
      particle_moments_.temperature.at(species)(elem_id) +=
        norm_squared * bias_corrected_weight * particle.species.mass /
        (velocity_dims_ * constants::boltzmann_constant * sum_of_weights_in_element);
    }
  }

  return this->particle_moments_.temperature;
}

std::unordered_map<Species, mfem::Vector>&
ParticleOperations::getVarianceReducedTemperature(
  const ParticleContainer& particles,
  const LowFidelityState* low_fidelity_state,
  const LowFidelityOperations* low_fidelity_operations
) {
  mfem::FiniteElementSpace& finite_element_space = discretization_.getFeSpace();
  mfem::Mesh& mesh = *finite_element_space.GetMesh();
  const int number_of_elements = finite_element_space.GetNE();

  const auto& low_fidelity_number_density_integral = control_variate_moment_integrals_.number_density;
  const auto& low_fidelity_bulk_velocity_integral = control_variate_moment_integrals_.momentum;
  const auto& low_fidelity_temperature_integral = control_variate_moment_integrals_.temperature;

  /*
   * Start with standard PIC temperature. Cells that successfully use
   * variance reduction will overwrite these values below.
   */
  for (auto& species_and_temperature : variance_reduced_particle_moments_.temperature) {
    const Species& species = species_and_temperature.first;
    species_and_temperature.second = particle_moments_.temperature.at(species);
  }

  /*
   * Uncorrected CV temperature and sum of squared residual weights.
   *
   * For particle i:
   *
   *   residual_weight_i = w_i * (1 - f_LF / f)
   */
  auto uncorrected_temperature = variance_reduced_particle_moments_.temperature;

  auto residual_weight_squared_sum = variance_reduced_particle_moments_.temperature;

  for (auto& species_and_temperature : uncorrected_temperature)
    species_and_temperature.second = 0.0;

  for (auto& species_and_sum : residual_weight_squared_sum)
    species_and_sum.second = 0.0;

  std::unordered_map<Species, std::vector<char>> variance_reduction_valid;

  variance_reduction_valid.reserve(variance_reduced_particle_moments_.temperature.size());

  for (const auto& species_and_temperature : variance_reduced_particle_moments_.temperature) {
    const Species& species = species_and_temperature.first;

    variance_reduction_valid.emplace(
      species,
      std::vector<char>(number_of_elements, 1)
    );
  }

  mfem::Vector fluctuation_velocity(velocity_dims_);

  for (const Particle& particle : particles) {
    if (!particle.is_alive)
      continue;

    const Species& species = particle.species;
    const int elem_id = particle.element;
    const double element_volume = mesh.GetElementVolume(elem_id);

    const mfem::Vector particle_velocity(particle.velocity.GetData(), velocity_dims_);

    if (!variance_reduction_elements_.at(species)[elem_id])
      continue;

    double low_fidelity_pdf = evaluateControlVariate_(particle, low_fidelity_state, low_fidelity_operations);
    const double high_fidelity_pdf = particle.particle_distribution_function_value;
    if (
      !(high_fidelity_pdf > 0.0) ||
      !std::isfinite(high_fidelity_pdf) ||
      !std::isfinite(low_fidelity_pdf)
    ) {
      variance_reduction_valid.at(species)[elem_id] = 0;
      continue;
    }

    const double residual_multiplier = 1.0 - low_fidelity_pdf / high_fidelity_pdf;
    const double residual_weight = particle.weight * residual_multiplier;

    if (!std::isfinite(residual_weight)) {
      variance_reduction_valid.at(species)[elem_id] = 0;
      continue;
    }

    double fluctuation_speed_squared = 0.0;

    for (int vel_dim = 0; vel_dim < velocity_dims_; ++vel_dim) {
      fluctuation_velocity(vel_dim) =
        particle_velocity(vel_dim) -
       variance_reduced_particle_moments_.bulk_velocity
          .at(species)(vel_dim, elem_id);

      fluctuation_speed_squared += fluctuation_velocity(vel_dim) * fluctuation_velocity(vel_dim);
    }

    const double corrected_number_density = variance_reduced_particle_moments_.number_density.at(species)(elem_id);
    const double corrected_mass = corrected_number_density * element_volume;

    if (
      !(corrected_mass > 0.0) ||
      !std::isfinite(corrected_mass)
    ) {
      variance_reduction_valid.at(species)[elem_id] = 0;
      continue;
    }

    const double m_over_dkb = species.mass / (velocity_dims_ * constants::boltzmann_constant);

    uncorrected_temperature.at(species)(elem_id) +=
      m_over_dkb *
      residual_weight *
      fluctuation_speed_squared /
      corrected_mass;

    residual_weight_squared_sum.at(species)(elem_id) +=
      residual_weight * residual_weight;
  }

  /*
   * Add the analytic LF contribution and then apply the cell-level
   * effective-weight correction to the entire central moment.
   */
  for (const auto& [species, elements] : variance_reduction_elements_) {
    for (int elem_id = 0; elem_id < number_of_elements; ++elem_id) {
      if (!elements[elem_id])
        continue;

      if (!variance_reduction_valid.at(species)[elem_id])
        continue;

      const double element_volume = mesh.GetElementVolume(elem_id);
      const double corrected_number_density = variance_reduced_particle_moments_.number_density.at(species)(elem_id);
      const double corrected_mass = corrected_number_density * element_volume;

      if (!(corrected_mass > 0.0) || !std::isfinite(corrected_mass)) {
        continue;
      }

      const double m_over_dkb = species.mass / (velocity_dims_ * constants::boltzmann_constant);
      double corrected_bulk_velocity_squared = 0.0;
      double corrected_velocity_dot_lf_momentum = 0.0;

      for (int vel_dim = 0; vel_dim < velocity_dims_; ++vel_dim) {
        const double corrected_velocity = variance_reduced_particle_moments_.bulk_velocity.at(species)(vel_dim, elem_id);
        corrected_bulk_velocity_squared += corrected_velocity * corrected_velocity;
        corrected_velocity_dot_lf_momentum += corrected_velocity * low_fidelity_bulk_velocity_integral.at(species)(vel_dim, elem_id);
      }

      double low_fidelity_temperature_contribution = low_fidelity_temperature_integral.at(species)(elem_id) / corrected_mass;
      low_fidelity_temperature_contribution += m_over_dkb * corrected_bulk_velocity_squared * low_fidelity_number_density_integral.at(species)(elem_id) / corrected_mass;
      low_fidelity_temperature_contribution -= 2.0 * m_over_dkb * corrected_velocity_dot_lf_momentum / corrected_mass;

      double& temperature = uncorrected_temperature.at(species)(elem_id);
      temperature += low_fidelity_temperature_contribution;
      const double squared_mass = corrected_mass * corrected_mass;
      const double squared_residual_weight_sum = residual_weight_squared_sum.at(species)(elem_id);
      const double bias_denominator = squared_mass - squared_residual_weight_sum;

      /*
       * A singular or negative denominator can occur because control
       * variate residual weights may be signed. In that case, leave
       * the output at its standard-PIC fallback value.
       */
      if (
        !(bias_denominator > 0.0) ||
        !std::isfinite(bias_denominator) ||
        !std::isfinite(temperature)
      ) {
        continue;
      }

      const double bias_factor = squared_mass / bias_denominator;
      const double corrected_temperature = bias_factor * temperature;
      if (!std::isfinite(corrected_temperature))
        continue;

      variance_reduced_particle_moments_.temperature.at(species)(elem_id) = corrected_temperature;
    }
  }

  return variance_reduced_particle_moments_.temperature;
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


bool ParticleOperations::hasControlVariate_(
  const Species& species,
  const LowFidelityState* low_fidelity_state
) const {
  if (variance_reduction_parameters_.strategy == VarianceReductionParameters::Strategy::UniformMaxwellian)
  {
    const std::vector<std::string>& reference_species = variance_reduction_parameters_.reference_species;
    return reference_species.empty() || std::ranges::find(reference_species, species.name) != reference_species.end();
  }
  return low_fidelity_state != nullptr && low_fidelity_state->getSpeciesIndex(species) >= 0;
}

double ParticleOperations::evaluateControlVariate_(
  const Particle& particle,
  const LowFidelityState* low_fidelity_state,
  const LowFidelityOperations* low_fidelity_operations
) const {
  const mfem::Vector particle_velocity(particle.velocity.GetData(), velocity_dims_);

  if (variance_reduction_parameters_.strategy == VarianceReductionParameters::Strategy::UniformMaxwellian)
  {
    const mfem::Vector reference_primitive_state{
      variance_reduction_parameters_.reference_number_density,
      variance_reduction_parameters_.reference_bulk_velocity[0],
      variance_reduction_parameters_.reference_bulk_velocity[1],
      variance_reduction_parameters_.reference_bulk_velocity[2],
      variance_reduction_parameters_.reference_temperature};
    return euler::evaluateMaxwellian(reference_primitive_state, particle_velocity, particle.species);
  }

  if (low_fidelity_state == nullptr || low_fidelity_operations == nullptr) {
    errorWithDeveloperMessage("Variance reduction without a low-fidelity model requires Strategy: Uniform Maxwellian.");
  }
  const mfem::Vector particle_position(particle.position.GetData(), dim_);
  return low_fidelity_operations->evaluateParticleDistributionFunction(
    *low_fidelity_state,
    particle_position,
    particle_velocity,
    particle.element,
    low_fidelity_state->getSpeciesIndex(particle.species));
}

std::unordered_map<Species, std::vector<char>> ParticleOperations::computeVarianceReductionElements_(
  const ParticleContainer& particles,
  const LowFidelityState* low_fidelity_state,
  const LowFidelityOperations* low_fidelity_operations
) const {
  const int number_of_elements = discretization_.getFeSpace().GetNE();

  std::unordered_map<Species, std::vector<char>> variance_reduction_elements;
  for (const Species& species : std::views::keys(sum_of_weights_)) {
    const char has_control_variate = hasControlVariate_(species, low_fidelity_state) ? 1 : 0;
    variance_reduction_elements.emplace(species, std::vector<char>(number_of_elements, has_control_variate));
  }

  if (!variance_reduction_parameters_.limit_variance_reduction)
    return variance_reduction_elements;

  for (const Particle& particle : particles) {
    if (not particle.is_alive) continue;

    char& use_variance_reduction = variance_reduction_elements.at(particle.species)[particle.element];
    if (!use_variance_reduction) continue;

    const double noise_reducing_factor =
      1.0 - evaluateControlVariate_(particle, low_fidelity_state, low_fidelity_operations) / particle.particle_distribution_function_value;
    // Written so that a NaN factor also disables variance reduction
    if (!(std::abs(noise_reducing_factor) < 1.0))
      use_variance_reduction = 0;
  }
  return variance_reduction_elements;
}

ParticleOperations::ControlVariateMomentIntegrals ParticleOperations::computeControlVariateMomentIntegrals_(
  const LowFidelityState* low_fidelity_state,
  const DGEulerOperations* low_fidelity_operations
) const {
  const mfem::FiniteElementSpace& finite_element_space = discretization_.getFeSpace();
  mfem::Mesh& mesh = *finite_element_space.GetMesh();
  const int number_of_elements = finite_element_space.GetNE();

  ControlVariateMomentIntegrals integrals;
  if (!(variance_reduction_parameters_.strategy == VarianceReductionParameters::Strategy::UniformMaxwellian))
  {
    if (low_fidelity_state == nullptr || low_fidelity_operations == nullptr)
      return integrals;
    integrals.number_density =
      low_fidelity_operations->integralForVarianceReducedNumberDensity(finite_element_space, *low_fidelity_state);
    integrals.momentum =
      low_fidelity_operations->integralForVarianceReducedBulkVelocity(finite_element_space, *low_fidelity_state);
    integrals.temperature =
      low_fidelity_operations->integralForVarianceReducedTemperature(finite_element_space, *low_fidelity_state, velocity_dims_);
    return integrals;
  }

  const double reference_number_density = variance_reduction_parameters_.reference_number_density;
  const std::vector<double>& reference_bulk_velocity = variance_reduction_parameters_.reference_bulk_velocity;
  double reference_bulk_velocity_squared = 0.0;
  for (int vel_dim = 0; vel_dim < velocity_dims_; ++vel_dim)
    reference_bulk_velocity_squared += reference_bulk_velocity[vel_dim] * reference_bulk_velocity[vel_dim];

  for (const Species& species : std::views::keys(sum_of_weights_)) {
    if (!hasControlVariate_(species, low_fidelity_state)) continue;

    const double m_over_dkb = species.mass / (velocity_dims_ * constants::boltzmann_constant);
    mfem::Vector& number_density = integrals.number_density.emplace(species, mfem::Vector(number_of_elements)).first->second;
    mfem::DenseMatrix& momentum = integrals.momentum.emplace(species, mfem::DenseMatrix(3, number_of_elements)).first->second;
    mfem::Vector& temperature = integrals.temperature.emplace(species, mfem::Vector(number_of_elements)).first->second;

    for (int elem_id = 0; elem_id < number_of_elements; ++elem_id) {
      const double element_number = reference_number_density * mesh.GetElementVolume(elem_id);
      number_density(elem_id) = element_number;
      for (int vel_dim = 0; vel_dim < 3; ++vel_dim)
        momentum(vel_dim, elem_id) = element_number * reference_bulk_velocity[vel_dim];
      temperature(elem_id) = element_number *
        (variance_reduction_parameters_.reference_temperature + m_over_dkb * reference_bulk_velocity_squared);
    }
  }
  return integrals;
}

mfem::Vector ParticleOperations::assembleControlVariateCharge_(
  const Species& species,
  const LowFidelityState* low_fidelity_state,
  const LowFidelityOperations* low_fidelity_operations
) const {
  if (variance_reduction_parameters_.strategy == VarianceReductionParameters::Strategy::UniformMaxwellian)
  {
    mfem::ConstantCoefficient charge_density(species.charge * variance_reduction_parameters_.reference_number_density);
    mfem::LinearForm control_variate_charge(&discretization_.getFeSpace());
    control_variate_charge.AddDomainIntegrator(new mfem::DomainLFIntegrator(charge_density));
    control_variate_charge.Assemble();
    return control_variate_charge;
  }

  return low_fidelity_operations->assembleChargePerSpecies(
    *low_fidelity_state,
    low_fidelity_state->getSpeciesIndex(species)).getIntegratedCharge();
}

} // namespace mfpic
