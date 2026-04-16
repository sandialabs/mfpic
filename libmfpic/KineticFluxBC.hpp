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
 * @brief A dummy class that will error out if any functions are called. 
 * This is needed because we are creating a mfem::NumericalFlux that is only 
 * used for boundary terms. mfem::NumericalFlux requires a FluxFunction even if we don't
 * use it.
 */

class DummyFluxFunction : public mfem::FluxFunction {
  
  public:

  DummyFluxFunction(const int num_equations, const int dim)
    : mfem::FluxFunction(num_equations, dim) {};

  double ComputeFlux(const mfem::Vector &, mfem::ElementTransformation &, mfem::DenseMatrix &) const override
  {
    errorWithDeveloperMessage("Evaluation routines should never be called!");
    return -1;
  }

};

/**
 * @brief A mfem::NumericalFlux that implements an absorbing (wall) 
 * boundary condition with a kinetic flux vector splitting approach.
 * See II.D in Phys. Plasmas 27, 113505 (2020)
 *
 * @warning Should (and can) only be used for boundary integration
 *
 * @note A mfem::NumericalFlux function is needed to pass to the
 * mfem::HyperbolicForm integrator. During boundary integration,
 * mfem will call Eval and we ignore the right data. 
 * \ref DummyFluxFunction should be used with this, 
 * which will throw an error if any functions are called. 
 *
 * @todo Assumes 5-moment fluid
 */

class KineticFluxNumericalFlux : public mfem::NumericalFlux {

  public:

  /**
   * @brief Construct a new KineticFlux with given spatial dimension for a \ref Species
   *
   * @param flux_function 
   * @param species Species type
   */

  KineticFluxNumericalFlux(const mfem::FluxFunction & flux_function, const Species species)
    : mfem::NumericalFlux(flux_function),
      species_(species) {}

  /**
   * @brief Compute normal flux
   *
   * @warning Only appropriate on boundary
   *
   * @todo Not sure if the characteristic speed makes sense
   *
   * @param conservative_state state at integration point
   * @param not_used right state is not meaningful
   * @param normal normal vector, usually not a unit vector
   * @param transformation element transformation with the integration point
   * @param flux_dot_n [out] storage for normal flux output 
   * @return maximum characteristic speed, fluid velocity plus speed of sound
   */
  double Eval(const mfem::Vector &conservative_state, const mfem::Vector &, const mfem::Vector &normal,
              mfem::FaceElementTransformations &transformation,
              mfem::Vector &flux_dot_n) const override;

  private:
    const Species species_;
};

/**
 * @brief A DG BC that uses the kinetic flux vector splitting approach
 */

struct KineticFluxBC : public DGBC {

  KineticFluxBC() = delete;
  KineticFluxBC(const int boundary_attribute, const mfem::Mesh& mesh, const Species species) 
  : DGBC(boundary_attribute, mesh, species) {};

  virtual ~KineticFluxBC() = default;

  std::unique_ptr<mfem::NonlinearFormIntegrator> makeIntegrator(const mfem::NumericalFlux & dg_assembly_numerical_flux) override {
    flux_function_ = std::make_unique<DummyFluxFunction>(euler::ConservativeVariables::NUM_VARS, dg_assembly_numerical_flux.GetFluxFunction().dim);
    numerical_flux_ = std::make_unique<KineticFluxNumericalFlux>(*flux_function_, species);

    return std::make_unique<mfem::HyperbolicFormIntegrator>(*numerical_flux_);
  }

  private:
    std::unique_ptr<DummyFluxFunction> flux_function_;
    std::unique_ptr<KineticFluxNumericalFlux> numerical_flux_;

};

} // namespace mfpic
