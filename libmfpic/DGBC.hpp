#pragma once

#include <libmfpic/Species.hpp>
#include <memory.h>
#include <mfem.hpp>
#include <mfem/fem/hyperbolic.hpp>

namespace mfpic {

 /**
 * @brief Parent struct for DG boundary conditions. 
 */
struct DGBC {

  DGBC() = delete;
  DGBC(const int boundary_attribute, const mfem::Mesh& mesh) 
  : boundary_attribute_has_boundary_condition(mesh.bdr_attributes.Max())
  {
    boundary_attribute_has_boundary_condition = false;
    boundary_attribute_has_boundary_condition[boundary_attribute - 1] = true;
  };

  virtual ~DGBC() = default;

  /**
   * @brief Generate a mfem::HyperbolicFormIntegrator that can be used as a BdrFaceIntegrator. This will apply the BC during assembly.
   *
   * @param dg_assembly_numerical_flux The numerical flux function from the assembly object this BC is being added to
   * @param species Species associated with the assembly object
   */
  virtual std::shared_ptr<mfem::HyperbolicFormIntegrator> makeIntegrator(const mfem::NumericalFlux & dg_assembly_numerical_flux, Species species) = 0;

  /// Boundary attribute numbers defining the boundaries to apply the BC
  mfem::Array<int> boundary_attribute_has_boundary_condition;

};

} // namespace
