#include <libmfpic/TextDataWriter.hpp>

#include <libmfpic/ElectrostaticFieldOperations.hpp>
#include <libmfpic/LowFidelityOperations.hpp>

#include <iomanip>

namespace mfpic {

namespace {

void writeHeaderToCSVFile(std::ofstream& csv_file) {
  csv_file << std::setprecision(std::numeric_limits<double>::max_digits10);
  csv_file << "# Time_Step Time Field_Energy Total_Fluid_Energy Total_Fluid_Kinetic_Energy Total_Charge" << std::endl;
}

void writeLineToCSVFile(
  std::ofstream& csv_file,
  const int i_time_step,
  const double time,
  const double field_energy,
  const double total_fluid_energy,
  const double total_fluid_kinetic_energy,
  const double total_charge)
{
  csv_file <<
  i_time_step << " " <<
  time << " " <<
  field_energy << " " <<
  total_fluid_energy << " " <<
  total_fluid_kinetic_energy << " " << 
  total_charge << std::endl;
}

}

TextDataWriter::TextDataWriter(const int num_low_fidelity_models,
                               const std::string &base_filename)
{
  std::string base = base_filename;
  const std::string ext = ".csv";
  if (base.size() >= ext.size() &&
      base.compare(base.size() - ext.size(), ext.size(), ext) == 0)
  {
    base.erase(base.size() - ext.size());
  }
  main_csv_file_.open(base + ".csv");
  writeHeaderToCSVFile(main_csv_file_);

  low_fidelity_csv_files_.clear();
  low_fidelity_csv_files_.resize(num_low_fidelity_models);
  for (int i = 0; i < num_low_fidelity_models; ++i)
  {
    low_fidelity_csv_files_[i].open(base + "_lf_" + std::to_string(i) + ".csv");
    writeHeaderToCSVFile(low_fidelity_csv_files_[i]);
  }
}



void TextDataWriter::output(
  const ElectrostaticFieldState& particle_electrostatic_field_state,
  const std::vector<ElectrostaticFieldState>& low_fidelity_field_states,
  const ElectrostaticFieldOperations& electrostatic_field_operations,
  const std::vector<LowFidelityState>& low_fidelity_states,
  const std::vector<std::unique_ptr<LowFidelityOperations>>& low_fidelity_operations_vector,
  const int i_time_step,
  const double time)
{
  const double particle_field_energy = electrostatic_field_operations.fieldEnergy(particle_electrostatic_field_state);
  const double particle_total_charge = electrostatic_field_operations.totalCharge(particle_electrostatic_field_state);
  writeLineToCSVFile(main_csv_file_, i_time_step, time, particle_field_energy, 0., 0.,particle_total_charge);

  for (int i_low_fidelity_model = 0; i_low_fidelity_model < std::ssize(low_fidelity_csv_files_); ++i_low_fidelity_model) {
    std::ofstream& low_fidelity_csv_file = low_fidelity_csv_files_[i_low_fidelity_model];

    const ElectrostaticFieldState& low_fidelity_field_state = low_fidelity_field_states[i_low_fidelity_model];
    const double low_fidelity_field_energy = electrostatic_field_operations.fieldEnergy(low_fidelity_field_state);

    const LowFidelityState& low_fidelity_state = low_fidelity_states[i_low_fidelity_model];
    const LowFidelityOperations& low_fidelity_operations = *low_fidelity_operations_vector[i_low_fidelity_model];
    const double total_fluid_energy = low_fidelity_operations.computeTotalEnergy(low_fidelity_state);
    const double total_fluid_kinetic_energy = low_fidelity_operations.computeTotalKineticEnergy(low_fidelity_state);
    const double total_fluid_charge = low_fidelity_operations.computeTotalCharge(low_fidelity_state);

    writeLineToCSVFile(
      low_fidelity_csv_file,
      i_time_step,
      time,
      low_fidelity_field_energy,
      total_fluid_energy,
      total_fluid_kinetic_energy,
      total_fluid_charge);
  }
}
}