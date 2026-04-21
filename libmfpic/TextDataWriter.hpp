#pragma once

#include <fstream>
#include <memory>
#include <vector>

namespace mfpic {

class ElectrostaticFieldOperations;
class ElectrostaticFieldState;
class LowFidelityOperations;
class LowFidelityState;

class TextDataWriter {
public:
  /**
   * @brief Construct a new Text Data Writer object
   */
  TextDataWriter(const int num_low_fidelity_models);

  /**
   * @brief Output data to the text files
   * 
   * @param electrostatic_field_state - the electrostatic potential
   * @param low_fidelity_field_states - the electrostatic potential from each low fidelity model being output
   * @param low_fidelity_states - the low fidelity states being output
   * @param i_time_step - the index of the timestep being output
   * @param time - the time for the output
   */
  void output(
    const ElectrostaticFieldState& electrostatic_field_state,
    const std::vector<ElectrostaticFieldState>& low_fidelity_field_states,
    const ElectrostaticFieldOperations& electrostatic_field_operations,
    const std::vector<LowFidelityState>& low_fidelity_states,
    const std::vector<std::unique_ptr<LowFidelityOperations>>& low_fidelity_operations,
    const int i_time_step,
    const double time);

private:
  std::ofstream main_csv_file_;
  std::vector<std::ofstream> low_fidelity_csv_files_;
};

}