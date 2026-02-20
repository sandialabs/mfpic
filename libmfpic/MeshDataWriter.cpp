#include <libmfpic/ElectrostaticFieldState.hpp>
#include <libmfpic/LowFidelityState.hpp>
#include <libmfpic/MeshDataWriter.hpp>


namespace mfpic {

MeshDataWriter::MeshDataWriter(const std::string& name, mfem::Mesh& mesh) : paraview_data_collection_(name, &mesh) {}

void MeshDataWriter::output(
  ElectrostaticFieldState& electrostatic_field_state,
  std::vector<LowFidelityState>& low_fidelity_states,
  const int i_time_step,
  const double time)
{
  paraview_data_collection_.SetCycle(i_time_step);
  paraview_data_collection_.SetTime(time);

  mfem::GridFunction potential_grid_function = electrostatic_field_state.getPotential();
  paraview_data_collection_.RegisterField("electrostatic_potential", &potential_grid_function);

  // TODO: the electric field doesn't appear quite right in ParaView, GetDerivative should be projecting onto an L2 finite element
  //  space not an HGRAD element space, also this is always being output to nodes in Paraview which is also skewing things.
  const int mesh_dim = paraview_data_collection_.GetMesh()->Dimension();

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

  if (low_fidelity_states.size() > 0) {
    LowFidelityState& low_fidelity_state = low_fidelity_states[0];
    for (int i_species = 0; i_species < low_fidelity_state.numSpecies(); ++i_species) {
      LowFidelitySpeciesState& species_state = low_fidelity_state.getSpeciesState(i_species);
      mfem::GridFunction& grid_function = species_state.getGridFunction();

      const std::string field_name = "species_" + std::to_string(i_species);
      paraview_data_collection_.RegisterField(field_name, &grid_function);
    }
  }

  paraview_data_collection_.Save();
}

}