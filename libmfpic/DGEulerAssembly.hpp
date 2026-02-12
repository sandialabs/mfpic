#pragma once

#include <libmfpic/DGAssembly.hpp>
#include <libmfpic/ElectromagneticFieldsEvaluator.hpp>
#include <libmfpic/EulerFlux.hpp>
#include <memory>
#include <mfem.hpp>

namespace mfpic {

class LowFidelitySpeciesState;
struct Species;

/** @brief Discontinous Galerkin operators for the Euler equations, which are of the form
 * \f[
 * \frac{\partial U}{\partial t} + \nabla \cdot F\left(U\right) = S\left(U\right)
 * \f]
 * where \f$F(U)\f$ is the hyperbolic flux function and \f$S(U)\f$ is the source term.
 *
 * @note the weak divergence operator is preassembled upon creation
 */

class DGEulerAssembly : public DGAssembly {

public:
 /**
  * @brief Construct a new DGEulerOperator.
  *
  * @param finite_element_space vector finite element space
  * @param species Species type
  */
  DGEulerAssembly(
    mfem::FiniteElementSpace &finite_element_space,
    const Species& species);
  
 /**
  * @brief Compute the electromagnetic source terms for the residual.
  * Adds \f$(S(U),v)\f$ to the vector \p rhs .
  *
  * @note Does not apply \f$M^{-1}\f$.
  *
  * @param species_state - current solution vector
  * @param field_evaluator electromagnetic field evaluator
  * @param[inout] rhs rhs storage
  */
  void computeSources(
    const LowFidelitySpeciesState& species_state,
    const ElectromagneticFieldsEvaluator& field_evaluator,
    mfem::Vector& rhs) const;

 /**
  * @brief Compute the integrated kinetic energy \f$(\mathcal{K},v)\f$ where
  * \f$\mathcal{K} = \frac{1}{2\rho} \bm{p} \cdot \bm{p}\f$ and add into the vector
  * the energy slot of the \p rhs .
  *
  * @note Does not apply \f$M^{-1}\f$ and the density and momentum slots are left unchanged.
  *
  * @param species_state current solution vector for a species
  * @param[inout] rhs rhs storage
  */
  void computeIntegratedKineticEnergy(const LowFidelitySpeciesState& species_state, mfem::Vector &rhs) const;

private:

  const mfem::real_t charge_over_mass_;

  DGEulerAssembly(
    mfem::FiniteElementSpace &finite_element_space,
    std::unique_ptr<EulerFlux> &&euler_flux,
    std::unique_ptr<mfem::RusanovFlux> &&numerical_euler_flux,
    mfem::real_t charge_over_mass) :
      DGAssembly(finite_element_space, std::move(euler_flux), std::move(numerical_euler_flux)),
      charge_over_mass_(charge_over_mass) {}

  static DGEulerAssembly CreateDGEulerAssembly_(
    mfem::FiniteElementSpace &finite_element_space,
    const Species & species)
  {
    auto euler_flux = std::make_unique<EulerFlux>(finite_element_space.GetMesh()->SpaceDimension(), species);
    auto numerical_euler_flux = std::make_unique<mfem::RusanovFlux>(*euler_flux);
    return DGEulerAssembly(finite_element_space, std::move(euler_flux), std::move(numerical_euler_flux), species.charge_over_mass);
  }
};
  
} // namespace mfpic
