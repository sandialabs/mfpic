#include <libmfpic/BGKRelaxationOperations.hpp>
#include <libmfpic/Constants.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/GenerateMaxwellianVelocity.hpp>
#include <libmfpic/ParticleOperations.hpp>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <random>
#include <unordered_set>

namespace mfpic {

BGKRelaxationOperations::BGKRelaxationOperations(double collision_frequency, const Species& species_to_relax,
  double averaging_timesteps, double smoothing_length, int smoothing_order,
  bool use_vr_moments, bool log_cell_energy, bool stamp_raw_density, bool ratio_stamp, int density_scale_stride, bool invariant_constraint, int kde_neighbors, bool relaxed_mass_renormalization, bool kde_unit_scale)
: collision_frequency_(collision_frequency),
  species_to_relax_(species_to_relax),
  averaging_timesteps_(averaging_timesteps),
  smoothing_length_(smoothing_length),
  smoothing_order_(smoothing_order),
  use_vr_moments_(use_vr_moments),
  log_cell_energy_(log_cell_energy),
  stamp_raw_density_(stamp_raw_density),
  ratio_stamp_(ratio_stamp),
  density_scale_stride_(density_scale_stride),
  invariant_constraint_(invariant_constraint),
  kde_neighbors_(kde_neighbors),
  relaxed_mass_renormalization_(relaxed_mass_renormalization),
  kde_unit_scale_(kde_unit_scale)
{
}

BGKRelaxationOperations::ConservedDensities BGKRelaxationOperations::smoothedConservedDensities_(
  const ParticleMoments& moments, const mfem::Mesh& const_mesh, int velocity_dims) const {
  mfem::Mesh& mesh = const_cast<mfem::Mesh&>(const_mesh);
  const mfem::Vector& number_densities = moments.number_density.at(species_to_relax_);
  const mfem::DenseMatrix& bulk_velocities = moments.bulk_velocity.at(species_to_relax_);
  const mfem::Vector& temperatures = moments.temperature.at(species_to_relax_);
  const mfem::Vector& bias_corrections = moments.temperature_bias_correction.at(species_to_relax_);
  const int num_elements = number_densities.Size();

  if (smoothing_neighbors_.empty()) {
    smoothing_neighbors_.resize(num_elements);
    std::vector<mfem::Vector> centers(num_elements);
    for (int element = 0; element < num_elements; ++element) mesh.GetElementCenter(element, centers[element]);
    for (int element = 0; element < num_elements; ++element)
      for (int other = 0; other < num_elements; ++other)
        if (centers[element].DistanceTo(centers[other]) <= smoothing_length_ * (1.0 + 1.0e-12))
          smoothing_neighbors_[element].push_back(other);
  }

  ConservedDensities local{mfem::Vector(num_elements), mfem::DenseMatrix(3, num_elements), mfem::Vector(num_elements)};
  local.number = 0.0; local.momentum = 0.0; local.speed_squared = 0.0;
  for (int element = 0; element < num_elements; ++element) {
    const double number_density = number_densities(element);
    if (not (number_density > 0.0)) continue;
    // Divide by the bias correction so n<|v|^2> is the particles' actual kinetic energy density.
    const double temperature = bias_corrections(element) > 0.0 ? temperatures(element) / bias_corrections(element) : 0.0;
    double bulk_speed_squared = 0.0;
    for (int dimension = 0; dimension < 3; ++dimension) {
      const double component = bulk_velocities(dimension, element);
      local.momentum(dimension, element) = number_density * component;
      if (dimension < velocity_dims) bulk_speed_squared += component * component;
    }
    local.number(element) = number_density;
    local.speed_squared(element) = number_density *
      (bulk_speed_squared + velocity_dims * constants::boltzmann_constant * (std::isfinite(temperature) ? temperature : 0.0) / species_to_relax_.mass);
  }

  const int space_dims = mesh.Dimension();
  const auto basis = [&](const mfem::Vector& offset) {
    std::vector<double> values{1.0};
    if (smoothing_order_ >= 1)
      for (int i = 0; i < space_dims; ++i) values.push_back(offset(i));
    if (smoothing_order_ >= 2)
      for (int i = 0; i < space_dims; ++i)
        for (int j = i; j < space_dims; ++j) values.push_back(offset(i) * offset(j));
    return values;
  };
  constexpr int num_quantities = 5;
  const auto quantity = [&](const ConservedDensities& densities, int q, int element) {
    if (q == 0) return densities.number(element);
    if (q == 4) return densities.speed_squared(element);
    return densities.momentum(q - 1, element);
  };
  const auto set_quantity = [&](ConservedDensities& densities, int q, int element, double value) {
    if (q == 0) densities.number(element) = value;
    else if (q == 4) densities.speed_squared(element) = value;
    else densities.momentum(q - 1, element) = value;
  };
  const auto valid = [&](const ConservedDensities& densities, int element) {
    const double n = densities.number(element);
    if (not (n > 0.0)) return false;
    double momentum_squared = 0.0;
    for (int dimension = 0; dimension < velocity_dims; ++dimension) momentum_squared += std::pow(densities.momentum(dimension, element), 2);
    return densities.speed_squared(element) * n - momentum_squared > 0.0;
  };

  ConservedDensities smoothed{mfem::Vector(num_elements), mfem::DenseMatrix(3, num_elements), mfem::Vector(num_elements)};
  smoothed.number = 0.0; smoothed.momentum = 0.0; smoothed.speed_squared = 0.0;
  ConservedDensities flat = smoothed;
  mfem::Vector center(space_dims), other_center(space_dims), offset(space_dims);
  for (int element = 0; element < num_elements; ++element) {
    mesh.GetElementCenter(element, center);
    const int num_basis = std::ssize(basis(offset));
    mfem::DenseMatrix normal_matrix(num_basis);
    normal_matrix = 0.0;
    std::vector<mfem::Vector> right_hand_sides(num_quantities, mfem::Vector(num_basis));
    for (mfem::Vector& rhs : right_hand_sides) rhs = 0.0;
    double total_volume = 0.0;
    for (const int other : smoothing_neighbors_[element]) {
      const double volume = mesh.GetElementVolume(other);
      total_volume += volume;
      mesh.GetElementCenter(other, other_center);
      subtract(other_center, center, offset);
      const std::vector<double> phi = basis(offset);
      for (int a = 0; a < num_basis; ++a) {
        for (int b = 0; b < num_basis; ++b) normal_matrix(a, b) += volume * phi[a] * phi[b];
        for (int q = 0; q < num_quantities; ++q) right_hand_sides[q](a) += volume * phi[a] * quantity(local, q, other);
      }
      for (int q = 0; q < num_quantities; ++q)
        set_quantity(flat, q, element, quantity(flat, q, element) + volume * quantity(local, q, other));
    }
    for (int q = 0; q < num_quantities; ++q) set_quantity(flat, q, element, quantity(flat, q, element) / total_volume);

    bool fit_ok = std::ssize(smoothing_neighbors_[element]) >= num_basis;
    if (fit_ok) {
      mfem::DenseMatrixInverse inverse(normal_matrix);
      mfem::Vector coefficients(num_basis);
      for (int q = 0; q < num_quantities; ++q) {
        inverse.Mult(right_hand_sides[q], coefficients);
        set_quantity(smoothed, q, element, coefficients(0));
      }
      fit_ok = valid(smoothed, element);
    }
    if (not fit_ok)
      for (int q = 0; q < num_quantities; ++q) set_quantity(smoothed, q, element, quantity(flat, q, element));
  }
  return smoothed;
}

void BGKRelaxationOperations::performCollisions(
  double dt,
  RandomNumberGenerator& generator,
  ParticleContainer& particles,
  ParticleOperations& particle_operations
) const {
  const ParticleMoments& particle_moments = particle_operations.getParticleMoments(particles);
  const int velocity_dims = particle_operations.getVelocityDims();
  const ConservedDensities current = smoothedConservedDensities_(particle_moments, particle_operations.getMesh(), velocity_dims);
  if (not averaged_densities_initialized_) {
    averaged_densities_ = current;
    averaged_densities_initialized_ = true;
  } else {
    const double alpha = 1.0 / averaging_timesteps_;
    averaged_densities_.number *= 1.0 - alpha;
    averaged_densities_.number.Add(alpha, current.number);
    averaged_densities_.momentum *= 1.0 - alpha;
    averaged_densities_.momentum.Add(alpha, current.momentum);
    averaged_densities_.speed_squared *= 1.0 - alpha;
    averaged_densities_.speed_squared.Add(alpha, current.speed_squared);
  }

  const int num_elements = current.number.Size();
  std::vector<mfem::Vector> primitive_states(num_elements);
  for (int element = 0; element < num_elements; ++element) {
    const double number_density = averaged_densities_.number(element);
    if (not (number_density > 0.0)) continue;
    mfem::Vector bulk_velocity(3);
    averaged_densities_.momentum.GetColumn(element, bulk_velocity);
    bulk_velocity /= number_density;
    double bulk_speed_squared = 0.0;
    for (int dimension = 0; dimension < velocity_dims; ++dimension) bulk_speed_squared += bulk_velocity(dimension) * bulk_velocity(dimension);
    const double temperature = species_to_relax_.mass * (averaged_densities_.speed_squared(element) / number_density - bulk_speed_squared)
      / (velocity_dims * constants::boltzmann_constant);
    primitive_states[element] = mfem::Vector{number_density, bulk_velocity(0), bulk_velocity(1), bulk_velocity(2), temperature};
  }
  if (density_scale_stride_ <= 0) {
    relaxToEquilibria_(dt, generator, particles, velocity_dims, primitive_states);
    return;
  }
  updateDensityScale_(averaged_densities_.number);
  relaxToEquilibria_(dt, generator, particles, velocity_dims, primitive_states, &density_scale_);
}

void BGKRelaxationOperations::updateDensityScale_(const mfem::Vector& number_densities) const {
  if (density_scale_.Size() != number_densities.Size()) {
    density_scale_ = number_densities;
    density_scale_sum_.SetSize(number_densities.Size());
    density_scale_sum_ = 0.0;
    density_scale_calls_ = 0;
  }
  density_scale_sum_ += number_densities;
  if (++density_scale_calls_ < density_scale_stride_) return;
  for (int element = 0; element < density_scale_.Size(); ++element)
    if (density_scale_sum_(element) > 0.0) density_scale_(element) = density_scale_sum_(element) / density_scale_calls_;
  density_scale_sum_ = 0.0;
  density_scale_calls_ = 0;
}

void BGKRelaxationOperations::relaxToEquilibria_(
  double dt,
  RandomNumberGenerator& generator,
  ParticleContainer& particles,
  int velocity_dims,
  const std::vector<mfem::Vector>& primitive_states,
  const mfem::Vector* stamp_number_densities
) const {
  const int num_elements = std::ssize(primitive_states);
  std::vector<std::vector<Particle*>> relaxing_particles_per_element =
    selectRelaxingParticlesPerElement_(dt, generator, particles, num_elements);

  cell_energy_change_.assign(num_elements, 0.0);
  for (int element = 0; element < num_elements; ++element) {
    const mfem::Vector& primitive_state = primitive_states[element];
    if (primitive_state.Size() != 5) continue;
    const double temperature = primitive_state(4);
    mfem::Vector bulk_velocity{primitive_state(1), primitive_state(2), primitive_state(3)};
    if (not (primitive_state(0) > 0.0 and std::isfinite(temperature) and temperature > 0.0 and bulk_velocity.CheckFinite() == 0)) continue;
    mfem::Vector stamp_state = primitive_state;
    if (stamp_number_densities != nullptr) stamp_state(0) = (*stamp_number_densities)(element);
    std::vector<double> equilibrium_masses_before;
    std::vector<mfem::Vector> velocities_before;
    if (invariant_constraint_)
      for (const Particle* particle : particles_per_element_[element]) {
        equilibrium_masses_before.push_back(particle->weight * euler::evaluateMaxwellian(stamp_state,
          mfem::Vector(particle->velocity.GetData(), velocity_dims), species_to_relax_) / particle->particle_distribution_function_value);
        velocities_before.push_back(particle->velocity);
      }
    const double thermal_speed = std::sqrt(constants::boltzmann_constant * temperature / species_to_relax_.mass);
    const auto normalized_velocity = [&](const mfem::Vector& velocity) {
      mfem::Vector xi(velocity_dims);
      for (int dimension = 0; dimension < velocity_dims; ++dimension) xi(dimension) = (velocity(dimension) - stamp_state(1 + dimension)) / thermal_speed;
      return xi;
    };
    std::vector<mfem::Vector> kde_velocities;
    std::vector<double> kde_local_weights, kde_particle_weights;
    if (kde_neighbors_ > 0)
      for (const Particle* particle : particles_per_element_[element]) {
        const mfem::Vector velocity(particle->velocity.GetData(), velocity_dims);
        const double local_weight = euler::evaluateMaxwellian(stamp_state, velocity, species_to_relax_) / particle->particle_distribution_function_value;
        if (not (std::isfinite(local_weight) and local_weight > 0.0)) continue;
        kde_velocities.push_back(normalized_velocity(velocity));
        kde_local_weights.push_back(local_weight);
        kde_particle_weights.push_back(particle->weight);
      }
    double relaxed_mass_before = 0.0, relaxed_mass_after = 0.0, relaxed_weight = 0.0;
    for (Particle* particle : relaxing_particles_per_element[element]) {
      const double energy_before = particle->velocity * particle->velocity;
      const double maxwellian_before = euler::evaluateMaxwellian(stamp_state,
        mfem::Vector(particle->velocity.GetData(), velocity_dims), species_to_relax_);
      relaxed_mass_before += particle->weight * maxwellian_before / particle->particle_distribution_function_value;
      particle->velocity = generateMaxwellianVelocity(bulk_velocity, temperature, species_to_relax_.mass, generator, velocity_dims);
      const mfem::Vector particle_velocity(particle->velocity.GetData(), velocity_dims);
      const double maxwellian_after = euler::evaluateMaxwellian(stamp_state, particle_velocity, species_to_relax_);
      const double ratio_stamped_value = particle->particle_distribution_function_value * maxwellian_after / maxwellian_before;
      const bool use_ratio = ratio_stamp_ and std::isfinite(ratio_stamped_value) and ratio_stamped_value > 0.0;
      particle->particle_distribution_function_value = use_ratio ? ratio_stamped_value : maxwellian_after;
      ratio_stamp_fallbacks_ += ratio_stamp_ and not use_ratio;
      if (kde_neighbors_ > 0 and not kde_local_weights.empty()) {
        const double local_weight = kdeLocalWeight_(normalized_velocity(particle_velocity), kde_velocities, kde_local_weights, kde_particle_weights);
        particle->particle_distribution_function_value = maxwellian_after / local_weight;
      }
      relaxed_mass_after += particle->weight * maxwellian_after / particle->particle_distribution_function_value;
      relaxed_weight += particle->weight;
      cell_energy_change_[element] += 0.5 * particle->weight * species_to_relax_.mass * (particle->velocity * particle->velocity - energy_before);
    }
    const double mass_ratio = relaxed_mass_after / (kde_unit_scale_ ? relaxed_weight : relaxed_mass_before);
    if ((relaxed_mass_renormalization_ or kde_unit_scale_) and std::isfinite(mass_ratio) and mass_ratio > 0.0)
      for (Particle* particle : relaxing_particles_per_element[element]) particle->particle_distribution_function_value *= mass_ratio;
    else if ((relaxed_mass_renormalization_ or kde_unit_scale_) and not relaxing_particles_per_element[element].empty())
      ++relaxed_mass_renormalization_failures_;
    if (invariant_constraint_)
      constrainInvariants_(element, stamp_state, velocity_dims, equilibrium_masses_before, velocities_before,
        relaxing_particles_per_element[element]);
  }
  if (invariant_constraint_)
    std::cout << "    BGK invariant constraint failures: " << invariant_constraint_failures_
              << " max Newton iterations: " << invariant_constraint_max_iterations_ << std::endl;
  if (relaxed_mass_renormalization_ or kde_unit_scale_)
    std::cout << "    BGK relaxed mass renormalization failures: " << relaxed_mass_renormalization_failures_ << std::endl;
  if (ratio_stamp_) std::cout << "    BGK ratio stamp fallbacks: " << ratio_stamp_fallbacks_ << std::endl;
}

void BGKRelaxationOperations::performCollisions(
  double dt,
  RandomNumberGenerator& generator,
  ParticleContainer& particles,
  ParticleOperations& particle_operations,
  const LowFidelityState& low_fidelity_state,
  const LowFidelityOperations& low_fidelity_operations
) const {
  if (not use_vr_moments_ or low_fidelity_state.getSpeciesIndex(species_to_relax_) < 0) {
    performCollisions(dt, generator, particles, particle_operations);
    return;
  }
  constexpr bool recompute_standard_pic_moments = true;
  const ParticleMoments& moments = particle_operations.getVarianceReducedParticleMoments(
    particles, low_fidelity_state, low_fidelity_operations, recompute_standard_pic_moments);
  const mfem::Vector& number_densities = moments.number_density.at(species_to_relax_);
  const mfem::DenseMatrix& bulk_velocities = moments.energy_conserving_bulk_velocity.at(species_to_relax_);
  const mfem::Vector& temperatures = moments.energy_conserving_temperature.at(species_to_relax_);
  std::vector<mfem::Vector> primitive_states(number_densities.Size());
  for (int element = 0; element < number_densities.Size(); ++element)
    primitive_states[element] = mfem::Vector{number_densities(element), bulk_velocities(0, element),
      bulk_velocities(1, element), bulk_velocities(2, element), temperatures(element)};

  std::vector<double> cell_kinetic_energy(number_densities.Size(), 0.0);
  if (log_cell_energy_)
    for (const Particle& particle : particles)
      if (particle.is_alive and particle.species == species_to_relax_)
        cell_kinetic_energy[particle.element] += 0.5 * particle.weight * species_to_relax_.mass * (particle.velocity * particle.velocity);

  mfem::Vector raw_number_densities;
  if (stamp_raw_density_) raw_number_densities = particle_operations.getParticleMoments(particles).number_density.at(species_to_relax_);
  relaxToEquilibria_(dt, generator, particles, particle_operations.getVelocityDims(), primitive_states,
    stamp_raw_density_ ? &raw_number_densities : nullptr);

  if (log_cell_energy_) {
    const std::vector<bool>& active_elements = particle_operations.getVarianceReductionElements().at(species_to_relax_);
    std::cout << std::setprecision(6) << "    BGK cell dE:";
    for (const double change : cell_energy_change_) std::cout << " " << change;
    std::cout << "\n    BGK cell E:";
    for (const double energy : cell_kinetic_energy) std::cout << " " << energy;
    std::cout << "\n    BGK cell active:";
    for (const bool active : active_elements) std::cout << " " << (active ? 1 : 0);
    std::cout << std::endl;
  }
}

double BGKRelaxationOperations::kdeLocalWeight_(
  const mfem::Vector& xi,
  const std::vector<mfem::Vector>& velocities,
  const std::vector<double>& local_weights,
  const std::vector<double>& particle_weights
) const {
  const int num_samples = std::ssize(velocities);
  std::vector<std::pair<double, int>> distances(num_samples);
  for (int j = 0; j < num_samples; ++j) distances[j] = {xi.DistanceSquaredTo(velocities[j]), j};
  const int num_neighbors = std::min(kde_neighbors_, num_samples);
  std::partial_sort(distances.begin(), distances.begin() + num_neighbors, distances.end());
  double weighted_sum = 0.0, weight_sum = 0.0;
  for (int n = 0; n < num_neighbors; ++n) {
    const int j = distances[n].second;
    weighted_sum += particle_weights[j] * local_weights[j];
    weight_sum += particle_weights[j];
  }
  return weighted_sum / weight_sum;
}

void BGKRelaxationOperations::constrainInvariants_(
  int element,
  const mfem::Vector& stamp_state,
  int velocity_dims,
  const std::vector<double>& equilibrium_masses_before,
  const std::vector<mfem::Vector>& velocities_before,
  const std::vector<Particle*>& relaxed_particles
) const {
  const std::vector<Particle*>& cell_particles = particles_per_element_[element];
  const std::unordered_set<const Particle*> relaxed_set(relaxed_particles.begin(), relaxed_particles.end());
  const int num_particles = std::ssize(cell_particles);
  const int num_constraints = 2 + velocity_dims;
  if (num_particles <= num_constraints) return;
  const double thermal_speed = std::sqrt(constants::boltzmann_constant * stamp_state(4) / species_to_relax_.mass);
  const auto invariants = [&](const mfem::Vector& velocity) {
    mfem::Vector values(num_constraints);
    values(0) = 1.0;
    values(num_constraints - 1) = 0.0;
    for (int dimension = 0; dimension < velocity_dims; ++dimension) {
      const double xi = (velocity(dimension) - stamp_state(1 + dimension)) / thermal_speed;
      values(1 + dimension) = xi;
      values(num_constraints - 1) += xi * xi;
    }
    return values;
  };

  mfem::Vector targets(num_constraints);
  targets = 0.0;
  for (int i = 0; i < num_particles; ++i) targets.Add(equilibrium_masses_before[i], invariants(velocities_before[i]));

  double relaxed_mass_before = 0.0, relaxed_weight = 0.0;
  for (int i = 0; i < num_particles; ++i)
    if (relaxed_set.contains(cell_particles[i])) {
      relaxed_mass_before += equilibrium_masses_before[i];
      relaxed_weight += cell_particles[i]->weight;
    }
  const double relaxed_mass_per_weight = relaxed_weight > 0.0 ? relaxed_mass_before / relaxed_weight : 0.0;

  std::vector<double> prior_masses(num_particles);
  std::vector<mfem::Vector> invariants_after(num_particles);
  for (int i = 0; i < num_particles; ++i) {
    const Particle& particle = *cell_particles[i];
    const double kde_mass = particle.weight * euler::evaluateMaxwellian(stamp_state,
      mfem::Vector(particle.velocity.GetData(), velocity_dims), species_to_relax_) / particle.particle_distribution_function_value;
    const bool relaxed = relaxed_set.contains(&particle);
    prior_masses[i] = not relaxed ? equilibrium_masses_before[i]
      : (kde_neighbors_ > 0 and std::isfinite(kde_mass) and kde_mass > 0.0) ? kde_mass
      : relaxed_mass_per_weight * particle.weight;
    invariants_after[i] = invariants(cell_particles[i]->velocity);
  }

  mfem::Vector multipliers(num_constraints), gradient(num_constraints), step(num_constraints);
  multipliers = 0.0;
  mfem::DenseMatrix hessian(num_constraints);
  std::vector<double> masses(num_particles);
  bool converged = false;
  int iteration = 0;
  for (; iteration < 50 and not converged; ++iteration) {
    gradient = targets;
    gradient.Neg();
    hessian = 0.0;
    for (int i = 0; i < num_particles; ++i) {
      masses[i] = prior_masses[i] * std::exp(multipliers * invariants_after[i]);
      gradient.Add(masses[i], invariants_after[i]);
      mfem::AddMult_a_VVt(masses[i], invariants_after[i], hessian);
    }
    double residual = 0.0;
    for (int k = 0; k < num_constraints; ++k) residual = std::max(residual, std::abs(gradient(k)));
    converged = residual <= 1.0e-12 * targets(0);
    if (converged or gradient.CheckFinite() != 0) break;
    mfem::DenseMatrixInverse inverse(hessian);
    inverse.Mult(gradient, step);
    multipliers -= step;
  }
  invariant_constraint_max_iterations_ = std::max(invariant_constraint_max_iterations_, iteration);
  const bool masses_valid = std::all_of(masses.begin(), masses.end(), [](double mass) { return std::isfinite(mass) and mass > 0.0; });
  if (not converged or multipliers.CheckFinite() != 0 or not masses_valid) {
    ++invariant_constraint_failures_;
    return;
  }
  for (int i = 0; i < num_particles; ++i) {
    Particle& particle = *cell_particles[i];
    particle.particle_distribution_function_value = particle.weight * euler::evaluateMaxwellian(stamp_state,
      mfem::Vector(particle.velocity.GetData(), velocity_dims), species_to_relax_) / masses[i];
  }
}

std::vector<std::vector<Particle*>> BGKRelaxationOperations::selectRelaxingParticlesPerElement_(
  double dt,
  RandomNumberGenerator& generator,
  ParticleContainer& particles,
  int num_elements
) const {
  std::uniform_real_distribution<double> uniform_unit_interval_distribution;
  const double relaxation_probability = -std::expm1(-collision_frequency_ * dt);

  std::vector<std::vector<Particle*>> relaxing_particles_per_element(num_elements);
  particles_per_element_.assign(num_elements, {});
  for (Particle& particle : particles) {
    if (not (particle.is_alive and particle.species == species_to_relax_)) continue;
    if (invariant_constraint_ or kde_neighbors_ > 0) particles_per_element_[particle.element].push_back(&particle);
    if (uniform_unit_interval_distribution(generator) < relaxation_probability) {
      relaxing_particles_per_element[particle.element].push_back(&particle);
    }
  }
  return relaxing_particles_per_element;
}

} // namespace mfpic
