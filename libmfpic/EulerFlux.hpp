#pragma once

#include <libmfpic/Euler.hpp>
#include <libmfpic/Species.hpp>
#include <mfem/mfem.hpp>

namespace mfpic {

/**
 * @brief A reimplementation of mfem::EulerFlux that is compatible
 * with a nD-3V paradigm
 */

class EulerFlux : public mfem::FluxFunction
{

public:
  /**
   * @brief Construct a new EulerFlux FluxFunction with given spatial
   * dimension for a \ref Species
   *
   * @param spatial_dim spatial dimension
   * @param species Species type
   */
  EulerFlux(const int spatial_dim, const Species species)
    : FluxFunction(euler::PrimitiveVariables::NUM_VARS, spatial_dim),
      species_(species) {}

  /**
   * @brief Compute Euler flux
   *
   * @param conservative_state state at current integration point
   * @param transformation current element transformation with the integration point
   * @param [out] flux storage for flux output 
   * @return maximum characteristic speed, fluid velocity plus speed of sound
   */
  double ComputeFlux(const mfem::Vector &conservative_state, mfem::ElementTransformation &transformation,
                     mfem::DenseMatrix &flux) const override;

  /**
   * @brief Compute normal flux
   *
   * @param conservative_state state at current integration point
   * @param normal normal vector, usually not a unit vector
   * @param transformation current element transformation with the integration point
   * @param flux_dot_n [out] storage for normal flux output 
   * @return maximum characteristic speed, fluid velocity plus speed of sound
   */
  double ComputeFluxDotN(const mfem::Vector &conservative_state, const mfem::Vector &normal,
                         mfem::FaceElementTransformations &transformation,
                         mfem::Vector &flux_dot_n) const override;

  private:

    const Species species_;

};

} // namespace mfpic
