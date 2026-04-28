#pragma once

#include <memory>

namespace YAML {
class Node;
}

namespace mfpic {

class Discretization;
class TimeIntegrator;

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

/**
 * @brief construct a time integrator based on user input
 * 
 * @param time_stepping_parameters - parameters defining time stepping
 * @param electrostatic_discretization - discretization of the electrostatic solve
 * @param push_low_fidelity_with_particle_fields - some time integrators switch which field state to use based on this flag
 * @return std::unique_ptr<TimeIntegrator> - pointer to the time integrator to be used in the simulation
 */
std::unique_ptr<TimeIntegrator> buildTimeIntegrator(
  const TimeSteppingParameters& time_stepping_parameters,
  Discretization& electrostatic_discretization,
  const bool push_low_fidelity_with_particle_fields);
}
