#include <libmfpic/ElectrostaticFieldState.hpp>
#include <libmfpic/LowFidelityState.hpp>
#include <libmfpic/MeshDataWriter.hpp>


namespace mfpic {

MeshDataWriter::MeshDataWriter(const std::string& name, mfem::Mesh& mesh) : paraview_data_collection_(name, &mesh) {}

void MeshDataWriter::output(
  ElectrostaticFieldState& electrostatic_field_state,
  std::vector<LowFidelityState>& /*low_fidelity_states*/,
  const int i_time_step,
  const double time)
{
  auto& message = std::cout;
  message << "MeshDataWrite::output" << std::endl;
  paraview_data_collection_.SetCycle(i_time_step);
  paraview_data_collection_.SetTime(time);

  message << "potential" << std::endl;
  mfem::GridFunction potential_grid_function = electrostatic_field_state.getPotential();
  paraview_data_collection_.RegisterField("electrostatic_potential", &potential_grid_function);

  // TODO: the electric field doesn't appear quite right in ParaView, GetDerivative should be projecting onto an L2 finite element
  //  space not an HGRAD element space, also this is always being output to nodes in Paraview which is also skewing things.
  const int mesh_dim = paraview_data_collection_.GetMesh()->Dimension();

  message << "electric field" << std::endl;
  std::vector<mfem::GridFunction> electric_field_grid_functions;
  electric_field_grid_functions.reserve(3);
  for (int i_dim = 0; i_dim < 3; ++i_dim) {
    electric_field_grid_functions.emplace_back(potential_grid_function.FESpace());
    if (i_dim < mesh_dim) {
      constexpr int component_to_differentiate = 1;
      potential_grid_function.GetDerivative(component_to_differentiate, i_dim, electric_field_grid_functions[i_dim]);
      electric_field_grid_functions[i_dim].Neg();
    } else {
      electric_field_grid_functions[i_dim] = 0;
    }
    paraview_data_collection_.RegisterField("E_" + std::to_string(i_dim), &electric_field_grid_functions[i_dim]);
  }

  // message << "low fidelity states" << std::endl;
  // if (low_fidelity_states.size() > 0) {
  //   message << "first low fidelity state" << std::endl;
  //   LowFidelityState& low_fidelity_state = low_fidelity_states[0];
  //   for (int i_species = 0; i_species < low_fidelity_state.numSpecies(); ++i_species) {
  //     message << "i_species = " << i_species << std::endl;
  //     LowFidelitySpeciesState& species_state = low_fidelity_state.getSpeciesState(i_species);
  //     mfem::GridFunction& grid_function = species_state.getGridFunction();
  //     message << "grid_function.Size() = " << grid_function.Size() << std::endl;

  //     message << "register field" << std::endl;
  //     const std::string field_name = "species_" + std::to_string(i_species);
  //     message << "&grid_function = " << &grid_function << std::endl;
  //     paraview_data_collection_.RegisterField(field_name, &grid_function);
  //   }
  // }

  message << "save" << std::endl;
  paraview_data_collection_.Save();
  message << "after save" << std::endl;
}

}