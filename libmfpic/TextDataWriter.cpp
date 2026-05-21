#include <libmfpic/TextDataWriter.hpp>

#include <libmfpic/ElectrostaticFieldOperations.hpp>
#include <libmfpic/LowFidelityOperations.hpp>

#include <iomanip>

namespace mfpic {

namespace {

void writeHeaderToCSVFile(std::ofstream& csv_file) {
  csv_file << std::setprecision(std::numeric_limits<double>::max_digits10);
  csv_file << "# Time_Step Time Field_Energy Total_Fluid_Energy Total_Fluid_Kinetic_Energy CFL" << std::endl;
}

void writeLineToCSVFile(
  std::ofstream& csv_file,
  const int i_time_step,
  const double time,
  const double field_energy,
  const double total_fluid_energy,
  const double total_fluid_kinetic_energy,
  const double cfl)
{
  csv_file <<
  i_time_step << " " <<
  time << " " <<
  field_energy << " " <<
  total_fluid_energy << " " <<
  total_fluid_kinetic_energy << " " <<
  cfl << std::endl;
}

}

TextDataWriter::TextDataWriter(const int num_low_fidelity_models) : main_csv_file_("output.csv") {
  writeHeaderToCSVFile(main_csv_file_);
  for (int i_low_fidelity_model = 0; i_low_fidelity_model < num_low_fidelity_models; ++i_low_fidelity_model) {
    low_fidelity_csv_files_.emplace_back("output_lf_" + std::to_string(i_low_fidelity_model) + ".csv");
    writeHeaderToCSVFile(low_fidelity_csv_files_[i_low_fidelity_model]);
  }
}

void TextDataWriter::output(
  const ElectrostaticFieldState& particle_electrostatic_field_state,
  const std::vector<ElectrostaticFieldState>& low_fidelity_field_states,
  const ElectrostaticFieldOperations& electrostatic_field_operations,
  const std::vector<LowFidelityState>& low_fidelity_states,
  const std::vector<std::unique_ptr<LowFidelityOperations>>& low_fidelity_operations_vector,
  const double timestep_size,
  const double smallest_cell_lengthscale,
  const int i_time_step,
  const double time)
{
  const double particle_field_energy = electrostatic_field_operations.fieldEnergy(particle_electrostatic_field_state);
  writeLineToCSVFile(main_csv_file_, i_time_step, time, particle_field_energy, 0., 0., 0.);

  for (int i_low_fidelity_model = 0; i_low_fidelity_model < std::ssize(low_fidelity_csv_files_); ++i_low_fidelity_model) {
    std::ofstream& low_fidelity_csv_file = low_fidelity_csv_files_[i_low_fidelity_model];

    const ElectrostaticFieldState& low_fidelity_field_state = low_fidelity_field_states[i_low_fidelity_model];
    const double low_fidelity_field_energy = electrostatic_field_operations.fieldEnergy(low_fidelity_field_state);

    const LowFidelityState& low_fidelity_state = low_fidelity_states[i_low_fidelity_model];
    const LowFidelityOperations& low_fidelity_operations = *low_fidelity_operations_vector[i_low_fidelity_model];
    const double total_fluid_energy = low_fidelity_operations.computeTotalEnergy(low_fidelity_state);
    const double total_fluid_kinetic_energy = low_fidelity_operations.computeTotalKineticEnergy(low_fidelity_state);
    const double cfl = low_fidelity_operations.estimateCFL(timestep_size, smallest_cell_lengthscale);

    writeLineToCSVFile(
      low_fidelity_csv_file,
      i_time_step,
      time,
      low_fidelity_field_energy,
      total_fluid_energy,
      total_fluid_kinetic_energy,
      cfl);
  }
}
}