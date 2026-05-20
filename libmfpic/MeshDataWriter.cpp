#include <libmfpic/ElectrostaticFieldState.hpp>
#include <libmfpic/LowFidelityState.hpp>
#include <libmfpic/MeshDataWriter.hpp>

namespace mfpic {

MeshDataWriter::MeshDataWriter(const std::string& name, mfem::Mesh& mesh) : paraview_data_collection_(name, &mesh) {}

void MeshDataWriter::output(
  ElectrostaticFieldState& particle_field_state,
  std::vector<ElectrostaticFieldState>& low_fidelity_field_states,
  std::vector<LowFidelityState>& low_fidelity_states,
  const int i_time_step,
  const double time)
{
  paraview_data_collection_.SetCycle(i_time_step);
  paraview_data_collection_.SetTime(time);

  const unsigned int num_lf_models = std::ssize(low_fidelity_states);

  std::vector<mfem::GridFunction> potential_grid_functions;
  potential_grid_functions.reserve(num_lf_models + 1);

  std::vector<mfem::GridFunction> electric_field_grid_functions;
  electric_field_grid_functions.reserve(3 * (num_lf_models + 1));

  std::vector<mfem::GridFunction> electric_field_dg_grid_functions;
  electric_field_dg_grid_functions.reserve(num_lf_models + 1);

  auto register_potential_and_e_field_for_model = [&](const std::string& suffix, ElectrostaticFieldState& field_state) {

    potential_grid_functions.emplace_back(field_state.getPotential());
    auto & potential_grid_function = potential_grid_functions.back();
    paraview_data_collection_.RegisterField("electrostatic_potential" + suffix, &potential_grid_function);

    electric_field_dg_grid_functions.emplace_back(field_state.getEFieldGridFunction());
    auto& e_field_dg_grid_function = electric_field_dg_grid_functions.back();
    paraview_data_collection_.RegisterField("e_field_dg" + suffix, &e_field_dg_grid_function);

    // TODO: the electric field doesn't appear quite right in ParaView, GetDerivative should be projecting onto an L2 finite element
    //  space not an HGRAD element space, also this is always being output to nodes in Paraview which is also skewing things.
    const int mesh_dim = paraview_data_collection_.GetMesh()->Dimension();

   for (int i_dim = 0; i_dim < 3; ++i_dim) {
      electric_field_grid_functions.emplace_back(potential_grid_function.FESpace());
      const int slot = electric_field_grid_functions.size() - 1;
      if (i_dim < mesh_dim) {
        constexpr int component_to_differentiate = 1;
        potential_grid_functions.back().GetDerivative(component_to_differentiate, i_dim, electric_field_grid_functions[slot]);
        electric_field_grid_functions[slot].Neg();
      } else {
        electric_field_grid_functions[slot] = 0;
      }
      paraview_data_collection_.RegisterField("E_" + std::to_string(i_dim) + suffix, &electric_field_grid_functions[slot]);
    }
  };

  register_potential_and_e_field_for_model("", particle_field_state);

  if (low_fidelity_states.size() > 0) {
    for (unsigned int i_lf_model = 0; i_lf_model < num_lf_models; ++i_lf_model) {
      LowFidelityState& low_fidelity_state = low_fidelity_states[i_lf_model];
      const std::string suffix = "_lf_" + std::to_string(i_lf_model);
      for (int i_species = 0; i_species < low_fidelity_state.numSpecies(); ++i_species) {
        LowFidelitySpeciesState& species_state = low_fidelity_state.getSpeciesState(i_species);
        mfem::GridFunction& grid_function = species_state.getGridFunction();

        const std::string field_name = "species_" + std::to_string(i_species) + suffix;
        paraview_data_collection_.RegisterField(field_name, &grid_function);
      }
      register_potential_and_e_field_for_model(suffix, low_fidelity_field_states[i_lf_model]);
    }
  }

  paraview_data_collection_.Save();
}

}
