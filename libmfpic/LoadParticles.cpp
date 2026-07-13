#include <libmfpic/LoadParticles.hpp>

namespace mfpic {

ParticleContainer loadParticles(
  const SourceParameters& source_parameters,
  RandomNumberGenerator& generator,
  std::shared_ptr<mfem::Mesh> mesh
) {
  ParticleContainer particles;

  auto number_density_function = [&] (const mfem::Vector& x) -> double {
    SourceStateParameters params_at_point = source_parameters.sourceStateParametersAtPoint(x);
    return params_at_point.number_density;
  };
  const mfem::Vector physical_particles_per_element = elementwiseIntegral(*mesh, number_density_function);
  const double total_physical_particles = std::reduce(
    physical_particles_per_element.begin(),
    physical_particles_per_element.end()
  );
  const double particle_weight = total_physical_particles / source_parameters.num_particles;

  MeshDistribution position_distribution(mesh, number_density_function);

  const Species& species = source_parameters.species;
  for (int i = 0; i < source_parameters.num_particles; ++i) {
    mfem::Vector position({0.0, 0.0, 0.0});
    const auto [random_mesh_position, element] = position_distribution.generateRandomPointAndElement(generator);
    for (int i = 0; i < random_mesh_position.Size(); i++) {
      position[i] = random_mesh_position[i];
    }

    SourceStateParameters source_state_parameters = source_parameters.sourceStateParametersAtPoint(random_mesh_position);

    // const double kappa = source_state_parameters.kappa;
    mfem::Vector primitive_state(5);
    primitive_state(0) = source_state_parameters.number_density;
    primitive_state(1) = source_state_parameters.bulk_velocity(0);
    primitive_state(2) = source_state_parameters.bulk_velocity(1);
    primitive_state(3) = source_state_parameters.bulk_velocity(2);
    primitive_state(4) = source_state_parameters.temperature;
    double particle_distribution_function_value = 0.0;

    mfem::Vector velocity = generate1DMaxwellianVelocity(
      source_state_parameters.bulk_velocity,
      source_state_parameters.temperature,
      species.mass,
      generator
    );
    particle_distribution_function_value = euler::evaluate1DMaxwellian(primitive_state,velocity,species);
    // if (kappa > 0.0) {
    //   velocity = generateIsotropicKappaVelocity(
    //     source_state_parameters.bulk_velocity,
    //     source_state_parameters.temperature,
    //     kappa,
    //     species.mass,
    //     generator
    //   );
    //   particle_distribution_function_value = euler::evaluateIsotropicKappaDistribution(primitive_state,velocity,kappa,species);
    // }
    // else {
    //   velocity = generateMaxwellianVelocity(
    //     source_state_parameters.bulk_velocity,
    //     source_state_parameters.temperature,
    //     species.mass,
    //     generator
    //   );
    //   particle_distribution_function_value = euler::evaluateMaxwellian(primitive_state,velocity,species);
    // }

    particles.addParticle(Particle{
      .position = position,
      .velocity = velocity,
      .element = element,
      .species = species,
      .weight = particle_weight,
      .is_alive = true,
      .particle_distribution_function_value = particle_distribution_function_value,
    });
  }

  return particles;
}

} // namespace mfpic
