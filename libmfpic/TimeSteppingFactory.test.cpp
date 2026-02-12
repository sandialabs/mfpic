#include <libmfpic/TimeSteppingFactory.hpp>

#include <yaml-cpp/yaml.h>

#include <gtest/gtest.h>


namespace {

using namespace mfpic;

TEST(TimeSteppingFactory, buildTimeSteppingParametersFromNumberOfTimestepsAndTimeStepSize) {
  constexpr int number_of_timesteps = 10;
  constexpr double timestep_size = 0.02;
  const std::string time_stepping_string(
    "Number of Time Steps: " + std::to_string(number_of_timesteps) + "\n"
    "Time Step Size: " + std::to_string(timestep_size) + "\n"
  );

  const YAML::Node time_stepping = YAML::Load(time_stepping_string);

  const TimeSteppingParameters time_stepping_parameters = buildTimeSteppingParametersFromYAML(time_stepping);

  EXPECT_EQ(number_of_timesteps, time_stepping_parameters.number_of_timesteps);
  EXPECT_EQ(timestep_size, time_stepping_parameters.timestep_size);
}

TEST(TimeSteppingFactory, buildTimeSteppingParametersFromNumberOfTimestepsAndFinalTime) {
  constexpr int number_of_timesteps = 5;
  constexpr double final_time = 0.02;
  const std::string time_stepping_string(
    "Number of Time Steps: " + std::to_string(number_of_timesteps) + "\n"
    "Final Time: " + std::to_string(final_time) + "\n"
  );

  const YAML::Node time_stepping = YAML::Load(time_stepping_string);

  const TimeSteppingParameters time_stepping_parameters = buildTimeSteppingParametersFromYAML(time_stepping);

  EXPECT_EQ(number_of_timesteps, time_stepping_parameters.number_of_timesteps);
  constexpr double timestep_size = final_time / number_of_timesteps;
  EXPECT_EQ(timestep_size, time_stepping_parameters.timestep_size);
}

TEST(TimeSteppingFactory, buildTimeSteppingParametersFromTimeStepSizeAndFinalTime) {
  constexpr double timestep_size = 0.011;
  constexpr double final_time = 0.2;
  const std::string time_stepping_string(
    "Time Step Size: " + std::to_string(timestep_size) + "\n"
    "Final Time: " + std::to_string(final_time) + "\n"
  );

  const YAML::Node time_stepping = YAML::Load(time_stepping_string);

  const TimeSteppingParameters time_stepping_parameters = buildTimeSteppingParametersFromYAML(time_stepping);

  EXPECT_LE(time_stepping_parameters.timestep_size, timestep_size);

  const double output_final_time = time_stepping_parameters.timestep_size * time_stepping_parameters.number_of_timesteps;
  EXPECT_DOUBLE_EQ(final_time, output_final_time);

  const double next_biggest_timestep_size = final_time / (time_stepping_parameters.number_of_timesteps - 1);
  EXPECT_LE(timestep_size, next_biggest_timestep_size);
}

}