#pragma once

#include <memory>
#include <unordered_map>
#include <vector>

namespace mfem {
  class Mesh;
}

namespace YAML {
  class Node;
}

namespace mfpic {

struct DGBC;
struct Species;

/**
 * @brief Parent struct for the parameters of an Euler fluid boundary condition.
 */
struct DGEulerBCParameters {
  virtual ~DGEulerBCParameters() = default;

  /**
   * @brief Build the boundary condition for one species on one boundary attribute.
   *
   * @param[in] boundary_attribute Boundary attribute on which to apply the BC.
   * @param[in] mesh Mesh containing the boundary.
   * @param[in] species Species of the Euler fluid.
   */
  virtual std::unique_ptr<DGBC> makeBC(int boundary_attribute, const mfem::Mesh& mesh, const Species& species) const = 0;
};

/**
 * @brief Parameters for a reflecting Euler fluid boundary condition.
 */
struct ReflectingDGEulerBCParameters : public DGEulerBCParameters {
  virtual std::unique_ptr<DGBC> makeBC(int boundary_attribute, const mfem::Mesh& mesh, const Species& species) const override;
};

/**
 * @brief Parameters for an absorbing Euler fluid boundary condition.
 */
struct AbsorbingDGEulerBCParameters : public DGEulerBCParameters {
  virtual std::unique_ptr<DGBC> makeBC(int boundary_attribute, const mfem::Mesh& mesh, const Species& species) const override;
};

/**
 * @brief Parameters for an Euler fluid boundary condition that absorbs outgoing fluid and
 * admits influx from a non-drifting Maxwellian immediately outside the boundary.
 */
struct MaxwellianInfluxDGEulerBCParameters : public DGEulerBCParameters {
  /**
   * @brief Ctor.
   *
   * @param[in] number_flux_density_in Flux density of incoming particles.
   * @param[in] temperature_in         Temperature of reservoir Maxwellian distribution.
   */
  MaxwellianInfluxDGEulerBCParameters(double number_flux_density_in, double temperature_in);

  virtual std::unique_ptr<DGBC> makeBC(int boundary_attribute, const mfem::Mesh& mesh, const Species& species) const override;

  /// Number flux density of particles from the external Maxwellian through boundary.
  const double number_flux_density;

  /// Temperature of external Maxwellian.
  const double temperature;
};

/**
 * @brief Parse the Euler fluid boundary conditions from the input deck
 *
 * @param[in] euler_fluids_bcs The "Boundary Conditions" sequence of the "Euler Fluids" block
 * @param[in] mesh_dimension Spatial dimension of the mesh, used to name inline mesh sides
 * @returns Map from boundary attribute to the parameters of the boundary condition on it
 */
std::unordered_map<int, std::unique_ptr<DGEulerBCParameters>> buildBoundaryAttributeToBCParametersFromYAML(
  const YAML::Node& euler_fluids_bcs,
  const int mesh_dimension);

/**
 * @brief Build boundary conditions for Euler system
 *
 * @returns Vector of vectors of \ref DGBC 's, where the outer vector is indexed by species consistent with \p species_list
 */
std::vector<std::vector<std::unique_ptr<DGBC>>> buildDGEulerBoundaryConditions(
  const std::unordered_map<int, std::unique_ptr<DGEulerBCParameters>>& boundary_attribute_to_bc_parameters,
  const mfem::Mesh& mesh,
  const std::vector<Species>& species_list);

}
