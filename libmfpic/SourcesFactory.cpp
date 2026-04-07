#include <libmfpic/Errors.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/SourcesFactory.hpp>

namespace mfpic {

namespace {
mfem::Vector constructPrimitiveState(const SourceStateParameters state) {
  return euler::constructPrimitiveState(state.number_density, state.bulk_velocity, state.temperature);
}
}

std::unique_ptr<mfem::VectorFunctionCoefficient> SourceParameters::getEulerVectorCoefficient() const {
  auto euler_vector_function = [this] (const mfem::Vector& x, mfem::Vector& conservative_state) {
    const SourceStateParameters source_state_parameters = sourceStateParametersAtPoint(x);
    const mfem::Vector primitive_state = constructPrimitiveState(source_state_parameters);
    conservative_state = euler::convertFromPrimitiveToConservative(primitive_state, species);
  };

  return std::make_unique<mfem::VectorFunctionCoefficient>(euler::ConservativeVariables::NUM_VARS, euler_vector_function);
}

ConstantSourceParameters::ConstantSourceParameters(
  const Species& species,
  const SourceStateParameters& state_parameters,
  const int num_particles)
  : SourceParameters(species, num_particles)
  , constant_state(state_parameters)
  {}

ConstantSourceParameters::ConstantSourceParameters(
  const Species& species,
  const double number_density,
  const double temperature,
  const mfem::Vector& bulk_velocity,
  const double kappa)
  : SourceParameters(species,kappa)
  , constant_state{.number_density = number_density, .bulk_velocity = bulk_velocity, .temperature = temperature, .kappa = kappa}
{}

SourceStateParameters ConstantSourceParameters::sourceStateParametersAtPoint(const mfem::Vector&) const {
  return constant_state;
}

SodSourceParameters::SodSourceParameters(
  const Species& species,
  const double discontinuity_location_in,
  const SourceStateParameters& left_state_parameters,
  const SourceStateParameters& right_state_parameters,
  const int num_particles)
  : SourceParameters(species, num_particles)
  , discontinuity_location(discontinuity_location_in)
  , left_state(left_state_parameters)
  , right_state(right_state_parameters)
{}

SourceStateParameters SodSourceParameters::sourceStateParametersAtPoint(const mfem::Vector& x) const {
  if (x[0] < discontinuity_location) {
    return left_state;
  } else {
    return right_state;
  }
}

GaussianSourceParameters::GaussianSourceParameters(
  const Species& species,
  const mfem::Vector& center,
  const double standard_deviation,
  const SourceStateParameters& offsets,
  const SourceStateParameters& heights,
  const double pressure_offset,
  const double pressure_height,
  const int num_particles)
  : SourceParameters(species, num_particles)
  , center(center)
  , standard_deviation(standard_deviation)
  , offsets(offsets)
  , heights(heights)
  , pressure_offset(pressure_offset)
  , pressure_height(pressure_height)
{}

SourceStateParameters GaussianSourceParameters::sourceStateParametersAtPoint(const mfem::Vector& x) const {
  mfem::Vector shifted_x(x.Size());
  for (int i_dim = 0; i_dim < x.Size(); ++i_dim) {
    shifted_x = x[i_dim] - center[i_dim];
  }
  const double exponential = exp(-0.5 * (shifted_x * shifted_x) / (standard_deviation * standard_deviation));

  const double number_density_offset = offsets.number_density;
  const mfem::Vector velocity_offset = offsets.bulk_velocity;
  const double number_density_height = heights.number_density;
  const mfem::Vector velocity_height = heights.bulk_velocity;

  const double number_density = number_density_offset + exponential * number_density_height;
  mfem::Vector velocity(velocity_offset);
  velocity.Add(exponential, velocity_height);
  const double pressure = pressure_offset + exponential * pressure_height;
  const double temperature = euler::temperature(number_density, pressure);

  return SourceStateParameters{
    .number_density = number_density,
    .bulk_velocity = velocity,
    .temperature = temperature,
  };
}

PeriodicPerturbationSourceParameters::PeriodicPerturbationSourceParameters(
  const Species& species,
  const mfem::Vector& wavevector,
  const SourceStateParameters& base_values,
  const SourceStateParameters& perturbations,
  const int num_particles)
  : SourceParameters(species, num_particles)
  , wavevector(wavevector)
  , base_values(base_values)
  , perturbations(perturbations)
{}

SourceStateParameters PeriodicPerturbationSourceParameters::sourceStateParametersAtPoint(const mfem::Vector& x) const {
  double k_dot_x = 0.;
  for (int i_dim = 0; i_dim < x.Size(); ++i_dim) {
    k_dot_x += wavevector[i_dim] * x[i_dim];
  }
  const double cos_kx = cos(k_dot_x);

  const double number_density = base_values.number_density * (1. + perturbations.number_density * cos_kx);
  mfem::Vector velocity(base_values.bulk_velocity.Size());
  for (int i_dim = 0; i_dim < velocity.Size(); ++ i_dim) {
    velocity[i_dim] = base_values.bulk_velocity[i_dim] * (1. + perturbations.bulk_velocity[i_dim] * cos_kx);
  }
  const double temperature = base_values.temperature * (1. + perturbations.temperature * cos_kx);

  return SourceStateParameters{
    .number_density = number_density,
    .bulk_velocity = velocity,
    .temperature = temperature,
  };
}

SourceStateParameters buildSourceStateParametersFromYAML(const YAML::Node& state_node) {
  const double number_density = state_node["Number Density"].as<double>();
  if (number_density <= 0.0) {
    errorWithUserMessage(formatParseMessage(state_node["Number Density"], "Number Density is nonpositive!"));
  }

  const YAML::Node& bulk_velocity_node = state_node["Bulk Velocity"];
  mfem::Vector bulk_velocity({0.0, 0.0, 0.0});
  if (bulk_velocity_node) {
    if (bulk_velocity_node.IsSequence()) {
      if (bulk_velocity_node.size() != 3) {
        errorWithUserMessage(formatParseMessage(bulk_velocity_node, "Bulk Velocity does not have length 3!"));
      }
      for (int i = 0; i < 3; i++) {
        bulk_velocity[i] = bulk_velocity_node[i].as<double>();
      }
    } else {
      std::cout << formatParseMessage(
        bulk_velocity_node, "Bulk Velocity is not a sequence, so it's ignored. Continuing.") << std::endl;
    }
  }

  double temperature = 0.;
  if (state_node["Temperature"]) {
    temperature = state_node["Temperature"].as<double>();
  }
  if (temperature < 0.0) {
    errorWithUserMessage(formatParseMessage(state_node["Temperature"], "Temperature is negative!"));
  }

  double kappa=-1;
  if (state_node["kappa"]) 
  {
    kappa = state_node["kappa"].as<double>();
    if (kappa <= 0.5) {
      errorWithUserMessage(formatParseMessage(state_node["kappa"], "kappa must be > 0.5!"));
    }
  }

  SourceStateParameters state_parameters{
    .number_density = number_density,
    .bulk_velocity = bulk_velocity,
    .temperature = temperature,
    .kappa = kappa
  };
  return state_parameters;
}

std::vector<std::unique_ptr<SourceParameters>> buildListOfSourceParametersFromYAML(
  const YAML::Node& sources_node,
  const std::unordered_map<std::string, Species>& species_map)
{
  assert(sources_node.IsSequence() or sources_node.IsNull() or not sources_node.IsDefined());

  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters;

  if (sources_node.IsSequence()) {
    for (const YAML::Node& source : sources_node) {
      std::vector<std::string> species_names;
      const YAML::Node& species_node = source["Species"];
      if (species_node.IsSequence()) {
        for (const YAML::Node& species : species_node) {
          const std::string species_name = species.as<std::string>();
          if (not species_map.contains(species_name)) {
            errorWithUserMessage(formatParseMessage(species, "Species was not created in the Species block!"));
          }
          species_names.push_back(species_name);
        }
      } else {
        errorWithUserMessage(formatParseMessage(species_node, "Species is not a sequence!"));
      }

      int num_particles_per_species = 0;
      if (source["Number of Macroparticles per Species"]) {
        num_particles_per_species = source["Number of Macroparticles per Species"].as<int>();
        if (num_particles_per_species <= 0) {
          errorWithUserMessage(formatParseMessage(
            source["Number of Macroparticles per Species"], "Number of Macroparticles per Species is nonpositive!"));
        }
      }

      if (source["Constant"]) {
        const YAML::Node& state_node = source["Constant"];
        const SourceStateParameters state_parameters = buildSourceStateParametersFromYAML(state_node);
        for (const std::string& species_name : species_names) {
          list_of_parameters.push_back(
            std::make_unique<ConstantSourceParameters>(
              species_map.at(species_name),
              state_parameters,
              num_particles_per_species));
        }
      } else if (source["Sod"]) {
        const YAML::Node& sod_node = source["Sod"];
        const double discontinuity_location = sod_node["Discontinuity Location"].as<double>();
        const YAML::Node& left_state_node = sod_node["Left State"];
        const SourceStateParameters left_state_parameters = buildSourceStateParametersFromYAML(left_state_node);

        const YAML::Node& right_state_node = sod_node["Right State"];
        const SourceStateParameters right_state_parameters = buildSourceStateParametersFromYAML(right_state_node);

        for (const std::string& species_name : species_names) {
          list_of_parameters.push_back(
            std::make_unique<SodSourceParameters>(
              species_map.at(species_name),
              discontinuity_location,
              left_state_parameters,
              right_state_parameters,
              num_particles_per_species));
        }
      } else if (source["Gaussian"]) {
        const YAML::Node& gaussian_node = source["Gaussian"];

        const YAML::Node& center_node = gaussian_node["Center"];
        mfem::Vector center(center_node.size());
        if (center_node.IsSequence()) {
          for (int i = 0; i < std::ssize(center_node); ++i){
            center[i] = center_node[i].as<double>();
          }
        } else {
          errorWithUserMessage(formatParseMessage(center_node, "Center must be a sequence.")); 
        }

        const double standard_deviation = gaussian_node["Standard Deviation"].as<double>();
        const SourceStateParameters offsets = buildSourceStateParametersFromYAML(gaussian_node["Offsets"]);
        const SourceStateParameters heights = buildSourceStateParametersFromYAML(gaussian_node["Heights"]);
        const double pressure_offset = gaussian_node["Offsets"]["Pressure"].as<double>();
        const double pressure_height = gaussian_node["Heights"]["Pressure"].as<double>();
        if (pressure_offset <= 0 or (pressure_height + pressure_offset) <= 0) {
          errorWithUserMessage(formatParseMessage(gaussian_node, "Pressure must be positive."));
        }

        for (const std::string& species_name : species_names) {
          list_of_parameters.push_back(std::make_unique<GaussianSourceParameters>(
            species_map.at(species_name),
            center,
            standard_deviation,
            offsets,
            heights,
            pressure_offset,
            pressure_height,
            num_particles_per_species));
        }
      } else if (source["Periodic Perturbation"]) {
        const YAML::Node& perturbation_node = source["Periodic Perturbation"];

        const SourceStateParameters base = buildSourceStateParametersFromYAML(perturbation_node["Base Values"]);
        const SourceStateParameters perturbations = buildSourceStateParametersFromYAML(perturbation_node["Perturbations"]);

        const YAML::Node& wavevector_node = perturbation_node["Wavevector"];
        mfem::Vector wavevector(wavevector_node.size());
        for (int i = 0; i < std::ssize(wavevector_node); ++i){
          wavevector[i] = wavevector_node[i].as<double>();
        }

        for (const std::string& species_name : species_names) {
          list_of_parameters.push_back(std::make_unique<PeriodicPerturbationSourceParameters>(
            species_map.at(species_name),
            wavevector,
            base,
            perturbations,
            num_particles_per_species));
        }
      } else {
        errorWithUserMessage(formatParseMessage(source, "It is required to either specify \"Constant\", \"Sod\", \"Gaussian\", or \"Periodic Perturbation\"."));
      }
    }
  } else {
    std::cout << formatParseMessage(
      sources_node, "Not a sequence, so it's ignored. Continuing.") << std::endl;
  }

  return list_of_parameters;
}

}
