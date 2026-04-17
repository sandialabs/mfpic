#include <libmfpic/TextDataWriter.hpp>

namespace mfpic {

namespace {

void writeHeaderToCSVFile(std::ofstream& csv_file) {
  csv_file << std::setprecision(std::numeric_limits<double>::digits);
  csv_file << "# Time_Step Time Field_Energy Total_Fluid_Energy" << std::endl;
}

void writeLineToCSVFile(std::ofstream& csv_file) {
  csv_file << i_time_step << " " << time << " " << field_energy << " " << total_fluid_energy << std::endl;
}

}

TextDataWriter::TextDataWriter(const int num_low_fidelity_models) : main_csv_file_("output.csv") {
  writeHeaderToCSVFile(main_csv_file_);
  for (int i_low_fidelity_model = 0; i_low_fidelity_model < num_low_fidelity_models; ++i_low_fidelity_model) {
    low_fidelity_csv_files_.emplace_back("output_lf_"+std::to_string(i)+".csv");
    writeHeaderToCSVFile(low_fidelity_csv_files_[i_low_fidelity_model]);
  }
}

void TextDataWriter::output(
  const ElectrostaticFieldState& particle_field_state,
  std::vector<ElectrostaticFieldState>& low_fidelity_field_states,
  std::vector<LowFidelityState>& low_fidelity_states,
  const int i_time_step,
  const double time)
{
  const particle_field_energy = electrostatic_field_operations->fieldEnergy(particle_electrostatic_field_state);
  writeLineToCSVFile(main_csv_file_, i_time_step, time, particle_field_energy, 0.);
  for (std::ofstream& low_fidelity_csv_file : low_fidelity_csv_files_) {
    writeLineToCSVFile(low_fidelity_csv_file, i_time_step, time, )
  }
}
}