#pragma once

namespace YAML {
class Node;
}

namespace mfpic {

enum class TimeIntegratorType { verlet, forward_euler };

struct TimeSteppingParameters {
  double timestep_size;
  int number_of_timesteps;

  TimeIntegratorType time_integrator_type = TimeIntegratorType::verlet;
};

/**
 * @brief Construct TimeSteppingParameters from yaml
 * 
 * @param time_stepping - yaml with 2 of the three keys "Number of Time Steps", "Time Step Size", and "Final Time"
 *  optionally specifies the time integrator type
 * @return TimeSteppingParameters 
 */
TimeSteppingParameters buildTimeSteppingParametersFromYAML(const YAML::Node& time_stepping);

}
