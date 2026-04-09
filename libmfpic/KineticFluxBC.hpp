#pragma once

#include <libmfpic/DGBC.hpp>
#include <libmfpic/Errors.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/Species.hpp>
#include <memory>
#include <mfem/fem/hyperbolic.hpp>
#include <mfem/mfem.hpp>

namespace mfpic {

/**
 * @brief A mfem::FluxFunction that implements an absorbing (wall) 
 * boundary condition with a kinetic flux vector splitting approach.
 * See II.D in Phys. Plasmas 27, 113505 (2020)
 *
 * @warning Should (and can) only be used for boundary integration
 *
 * @note A mfem::NumericalFlux function is needed to pass to the
 * mfem::HyperbolicForm integrator. During boundary integration,
 * by default mfem will only call FluxFunction::ComputeFluxDotN 
 * via the NumericalFlux. \ref KineticFluxBCNumericalFlux should be used, which
 * will throw an error if other functions are called. 
 *
 * @todo Assumes 5-moment fluid
 */

class KineticFluxFluxFunction : public mfem::FluxFunction {

  public:

  /**
   * @brief Construct a new KineticFluxBC with given spatial
   * dimension for a \ref Species
   *
   * @param spatial_dim spatial dimension
   * @param species Species type
   */

  KineticFluxFluxFunction(const int spatial_dim, const Species species)
    : mfem::FluxFunction(euler::PrimitiveVariables::NUM_VARS, spatial_dim),
      species_(species) {}
  /// Errors out if called
  double ComputeFlux (const mfem::Vector &, mfem::ElementTransformation &, mfem::DenseMatrix &) const override 
  { 
    errorWithDeveloperMessage("ComputeFlux cannot be called!");
    return -1;
  };

  /**
   * @brief Compute normal flux
   *
   * @warning Only appropriate on boundary
   *
   * @todo Not sure if the characteristic speed makes sense
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

/**
 * @brief A dummy class that will error out if anything other than
 * GetFluxFunction() is called. This is needed because mfem::HyperbolicIntegrator
 * requires a mfem::NumericalFlux and we only want to use the FluxFunction for
 * boundary assembly. The mfem::NumericalFlux is the only way to supply the FluxFunction
 * in the assembly routines.
 */

class DummyNumericalFlux : public mfem::NumericalFlux {
  
  public:

  DummyNumericalFlux(const mfem::FluxFunction & flux)
    : mfem::NumericalFlux(flux) {};

  double Eval(const mfem::Vector &, const mfem::Vector &, const mfem::Vector &,
              mfem::FaceElementTransformations &, mfem::Vector &) const override
  {
    errorWithDeveloperMessage("Evaluation routines should never be called!");
    return -1;
  }

};

struct KineticFluxBC : public DGBC {

  KineticFluxBC() = delete;
  KineticFluxBC(const int boundary_attribute, const mfem::Mesh& mesh) 
  : DGBC(boundary_attribute, mesh) {};

  virtual ~KineticFluxBC() = default;

  std::shared_ptr<mfem::HyperbolicFormIntegrator> makeIntegrator(const mfem::NumericalFlux & dg_assembly_numerical_flux, Species species) override {
    flux_function_ = std::make_unique<KineticFluxFluxFunction>(dg_assembly_numerical_flux.GetFluxFunction().dim, species);
    numerical_flux_ = std::make_unique<DummyNumericalFlux>(*flux_function_);

    return std::make_shared<mfem::HyperbolicFormIntegrator>(*numerical_flux_);
  }

  private:
    std::unique_ptr<KineticFluxFluxFunction> flux_function_;
    std::unique_ptr<DummyNumericalFlux> numerical_flux_;

};

} // namespace mfpic
