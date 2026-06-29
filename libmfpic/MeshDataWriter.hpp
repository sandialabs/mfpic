#pragma once

#include <mfem/mfem.hpp>

#include <string>

namespace mfpic {

class ElectrostaticFieldOperations;
class ElectrostaticFieldState;
class Discretization;
class LowFidelityOperations;
class LowFidelityState;

class MeshDataWriter {
public:
  /**
   * @brief Construct a new Mesh Data Writer object
   * 
   * @param folder_name - the name of the folder to contain the mesh data
   * @param mesh - the mesh that the data exists on, this must be nonconst because mfem::ParaViewDataCollection requires a
   *  nonconst mesh, but the mesh shouldn't be meaningfully changed by this object.
   */
  MeshDataWriter(const std::string& folder_name, mfem::Mesh& mesh);

  /**
   * @brief Output the electrostatic state and low fidelity state on the mesh
   * 
   * @param particle_field_state - the electrostatic potential being output
   * @param low_fidelity_field_states - the electrostatic potential from each low fidelity model being output
   * @param low_fidelity_states - the low fidelity states being output
   * @param i_time_step - the index of the timestep being output
   * @param time - the time for the output
   */
  void output(
    ElectrostaticFieldState& electrostatic_field_state,
    std::vector<ElectrostaticFieldState>& low_fidelity_field_states,
    ElectrostaticFieldOperations& electrostatic_field_operations,
    Discretization& electrostatic_discretization,
    std::vector<LowFidelityState>& low_fidelity_states,
    const std::vector<std::unique_ptr<LowFidelityOperations>>& low_fidelity_operations_vector,
    const int i_time_step,
    const double time);

private:
  mfem::ParaViewDataCollection paraview_data_collection_;
};

}
