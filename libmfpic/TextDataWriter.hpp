#pragma once

namespace mfpic {

class ElectrostaticFieldState;
class LowFidelityState;

class TextDataWriter {
public:
  /**
   * @brief Construct a new Text Data Writer object
   */
  TextDataWriter();

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
    ElectrostaticFieldState& electrostatic_field_state,
    std::vector<ElectrostaticFieldState>& low_fidelity_field_states,
    std::vector<LowFidelityState>& low_fidelity_states,
    const int i_time_step,
    const double time);

private:
  std::ofstream main_csv_file_;
  std::vector<std::ofstream> low_fidelity_csv_files_;
};

}