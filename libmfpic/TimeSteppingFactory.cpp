#include <libmfpic/TimeSteppingFactory.hpp>
#include <libmfpic/Errors.hpp>

#include <yaml-cpp/yaml.h>

namespace mfpic {

TimeSteppingParameters buildTimeSteppingParametersFromYAML(const YAML::Node& time_stepping) {
  TimeSteppingParameters parameters;

  if (time_stepping["Number of Time Steps"]) {
    const int number_of_timesteps = time_stepping["Number of Time Steps"].as<int>();
    parameters.number_of_timesteps = number_of_timesteps;
    if (time_stepping["Time Step Size"]) {
      const double timestep_size = time_stepping["Time Step Size"].as<double>();
      parameters.timestep_size = timestep_size;
    } else if (time_stepping["Final Time"]) {
      const double final_time = time_stepping["Final Time"].as<double>();
      const double timestep_size = final_time / number_of_timesteps;
      parameters.timestep_size = timestep_size;
    } else {
      errorWithUserMessage("Time Stepping must specify either \"Time Step Size\" or \"Final Time\"");
    }
  } else {
    const double timestep_size = time_stepping["Time Step Size"].as<double>();
    const double final_time = time_stepping["Final Time"].as<double>();

    const int number_of_timesteps = std::ceil(final_time / timestep_size);
    const double new_timestep_size = final_time / number_of_timesteps;
    parameters.number_of_timesteps = number_of_timesteps;
    parameters.timestep_size = new_timestep_size;
  }

  parameters.time_integrator_type = TimeIntegratorType::verlet;
  if (time_stepping["Type"]) {
    const std::string type_string = time_stepping["Type"].as<std::string>();
    if (type_string == "Verlet") {
      parameters.time_integrator_type = TimeIntegratorType::verlet;
    } else if (type_string == "Forward Euler") {
      parameters.time_integrator_type = TimeIntegratorType::forward_euler;
    } else {
      errorWithUserMessage("Time Stepping: Type: must be either \"Verlet\" or \"Forward Euler\".");
    }
  }

  return parameters;
}

} // namespace mfpic
