#pragma once

#include <libmfpic/KineticFluxBC.hpp>
#include <libmfpic/Species.hpp>

#include <mfem/mfem.hpp>

#include <memory>

namespace mfpic {

/**
 * @brief A \ref KineticFluxNumericalFlux that also adds the inward half-range flux of a
 * non-drifting Maxwellian immediately outside the boundary.
 *
 * @warning Should (and can) only be used for boundary integration.
 */
class MaxwellianInfluxNumericalFlux : public KineticFluxNumericalFlux {
public:
  /**
   * @brief Construct a new MaxwellianInfluxNumericalFlux for a \ref Species
   *
   * @param[in] flux_function Dummy flux function that the numerical flux must reference.
   * @param[in] species Species.
   * @param[in] number_flux_density Number flux density of particles from the external Maxwellian through boundary.
   * @param[in] temperature Temperature of external Maxwellian.
   */
  MaxwellianInfluxNumericalFlux(
    const mfem::FluxFunction& flux_function,
    const Species species,
    double number_flux_density,
    double temperature
  );

  /**
   * @brief Compute normal flux: outgoing kinetic flux of the interior state plus the external Maxwellian influx.
   *
   * @warning Only appropriate on boundary.
   *
   * @param[in]     conservative_state State at integration point.
   * @param[in]     not_used           "Right" state is not meaningful on boundaries.
   * @param[in]     normal             Outward normal vector, usually not a unit vector.
   * @param[in,out] transformation     Element transformation with the integration point.
   * @param[out]    flux_dot_n         Storage for normal flux output.
   *
   * @return Maximum characteristic speed of the interior state.
   */
  double Eval(
    const mfem::Vector &conservative_state,
    const mfem::Vector &not_used,
    const mfem::Vector &normal,
    mfem::FaceElementTransformations &transformation,
    mfem::Vector &flux_dot_n
  ) const override;

private:
  /// Number flux density of particles from the external Maxwellian through boundary.
  const double number_flux_density_;

  /// Temperature of external Maxwellian.
  const double temperature_;
};

/**
 * @brief A DG BC that absorbs outgoing fluid with the kinetic flux vector splitting approach and
 * admits influx from a non-drifting Maxwellian immediately outside the boundary
 */
struct MaxwellianInfluxBC : public KineticFluxBC {
  MaxwellianInfluxBC() = delete;

  /**
   * @brief Ctor.
   *
   * @param[in] boundary_attribute     MFEM attribute of boundary to apply condition.
   * @param[in] mesh                   Mesh.
   * @param[in] species                Species to which this boundary condition applies.
   * @param[in] number_flux_density_in Number flux density of physical particles entering through boundary.
   * @param[in] temperature_in         Temperature of Maxwellian immediately opposite boundary.
   */
  MaxwellianInfluxBC(
    int boundary_attribute,
    const mfem::Mesh& mesh,
    const Species& species,
    double number_flux_density_in,
    double temperature_in
  );

  /// Number flux density of particles from the external Maxwellian through boundary.
  const double number_flux_density;

  /// Temperature of external Maxwellian.
  const double temperature;

protected:
  /**
   * @brief Make the numerical flux applied on the boundary.
   *
   * @param flux_function Dummy flux function that the numerical flux must reference.
   */
  virtual std::unique_ptr<KineticFluxNumericalFlux> makeNumericalFlux_(const mfem::FluxFunction& flux_function) const override;
};

} // namespace mfpic
