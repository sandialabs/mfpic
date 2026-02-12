#pragma once

namespace YAML {
class Node;
}

namespace mfpic {

struct TimeSteppingParameters {
  double timestep_size;
  int number_of_timesteps;
};

/**
 * @brief Construct TimeSteppingParameters from yaml
 * 
 * @param time_stepping - yaml with 2 of the three keys "Number of Time Steps", "Time Step Size", and "Final Time"
 * @return TimeSteppingParameters 
 */
TimeSteppingParameters buildTimeSteppingParametersFromYAML(const YAML::Node& time_stepping);

}
