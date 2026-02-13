#pragma once

#include <libmfpic/DGGhostBoundaryIntegrator.hpp>
#include <libmfpic/ElectromagneticFieldsEvaluator.hpp>
#include <memory>
#include <mfem.hpp>
#include <mfem/config/config.hpp>
#include <mfem/fem/eltrans.hpp>
#include <mfem/fem/fe/fe_base.hpp>
#include <mfem/fem/fespace.hpp>
#include <mfem/fem/gridfunc.hpp>
#include <mfem/fem/hyperbolic.hpp>
#include <mfem/fem/intrules.hpp>
#include <mfem/fem/linearform.hpp>
#include <mfem/fem/lininteg.hpp>
#include <mfem/linalg/densemat.hpp>

namespace mfpic {

class LowFidelitySpeciesState;

/** @brief General Discontinous Galerkin assembly routines for solving hyperbolic equations of the form 
 * \f[
 * \frac{\partial U}{\partial t} + \nabla \cdot F\left(U\right) = S\left(U\right)
 * \f]
 * where \f$F(U)\f$ is the hyperbolic flux function and \f$S(U)\f$ is the source term.
 *
 */
class DGAssembly{

public:
 /**
  * @brief Construct a new DGAssembly.
  *
  * @param finite_element_space vector finite element space
  * @param flux_function flux function 
  * @param numerical_flux_function numerical flux function
  */
  DGAssembly(
    mfem::FiniteElementSpace &finite_element_space,
    std::unique_ptr<mfem::FluxFunction> &&flux_function,
    std::unique_ptr<mfem::NumericalFlux> &&numerical_flux_function
  );
  
  /**
   * @brief Compute the hyperbolic flux vector in an element
   *
   * @param element element index 
   * @param state_in_current_element solution vector in the element
   * @param[out] flux_in_current_element flux vector \f$F(U)\f$ in the element
   */
  void computeFluxVectorInElement(const int element, const mfem::DenseMatrix &state_in_current_element, mfem::DenseMatrix &flux_in_current_element) const;

  /**
   * @brief Applies the weak divergence operator to the hyperbolic flux vector and adds the results to the right hand side.
   *        Adds \f$\left(F(U), \nabla v\right)\f$ to \p rhs .
   *
   * @param element element index 
   * @param flux_in_current_element flux vector \f$F(U)\f$ in the element
   * @param[inout] rhs rhs storage 
   */
  void applyWeakDivergenceInElement(const int element, const mfem::DenseMatrix &flux_in_current_element, mfem::DenseMatrix &rhs_in_current_element) const;

  /**
   * @brief Compute the hyperbolic flux terms on the element boundaries, \f$\left\langle\hat{F} \cdot n, v\right\rangle\f$ 
   * and adds to \p rhs .
   *
   * @param dofs current solution vector
   * @param[inout] rhs rhs storage 
   */
  void computeFluxOnElementBoundaries(const mfem::Vector &dofs, mfem::Vector &rhs) const;

  /**
   * @brief Compute the hyperbolic flux terms for the residual.
   * Adds \f$\left(F(U), \nabla v\right) - \left\langle\hat{F} \cdot n, v\right\rangle\f$ to \p rhs .
   *
   * @note Does not apply \f$M^{-1}\f$.
   *
   * @param dofs current solution vector
   * @param[inout] rhs rhs storage
   */
  void computeHyperbolicFluxes(const mfem::Vector &dofs, mfem::Vector &rhs) const;

 /**
  * @brief Compute the electromagnetic source terms for the residual.
  * Adds \f$(S(U),v)\f$ to the vector \p rhs .
  *
  * @note Does not apply \f$M^{-1}\f$.
  *
  * @param species_state - current solution for a given species
  * @param field_evaluator electromagnetic field evaluator
  * @param[inout] rhs rhs storage
  */
  virtual void computeSources(
    const LowFidelitySpeciesState& species_state,
    const ElectromagneticFieldsEvaluator& field_evaluator,
    mfem::Vector& rhs) const = 0;

  /**
  * @brief Apply the local inverse mass matrices and store in \p dofs
  *
  * @param values Function values
  * @param[out] dofs Coefficient values
  */
  void applyInverseMass(const mfem::Vector &values, mfem::Vector &dofs) const;

  /**
   * @brief Make a new mfem::Vector sized for this operator's finite element space.
   */
  mfem::Vector getFEVector() const { return mfem::Vector(finite_element_space_.GetTrueVSize()); };

  mfem::real_t getNumberOfEquations() const { return num_equations_;}

  mfem::FiniteElementSpace & getFiniteElementSpace() const { return finite_element_space_;}

  /**
   * @brief Add a ghost cell boundary condition to the nonlinear form
   *
   * @param boundary_condition \ref DGGhostBC that sets the ghost cell state
   */
  void addGhostBoundaryCondition(std::unique_ptr<DGGhostBC> && boundary_condition)
  {
    nonlinear_form_->UseExternalIntegrators();
    ghost_bcs_.push_back(std::move(boundary_condition));
    const auto & current_bc = ghost_bcs_.back();
    ghost_bc_integrators_.push_back(
      std::make_unique<DGGhostBoundaryIntegrator>(*numerical_flux_function_, *current_bc));
    nonlinear_form_->AddBdrFaceIntegrator(
      (ghost_bc_integrators_.back()).get(), current_bc->boundary_attribute_has_boundary_condition);
  }

  /// get global maximum characteristic speed to be used in CFL condition
  /// where max_char_speed is updated during RHS evaluation
  mfem::real_t getMaxCharSpeed() const { return max_characteristic_speed_; }

  /// Dtor.
  virtual ~DGAssembly();

private:
  /// Spatial dimension
  const int dim_;

  /// Finite element space
  mfem::FiniteElementSpace &finite_element_space_;

  /// Base Nonlinear Form
  std::shared_ptr<mfem::NonlinearForm> nonlinear_form_;
  /// element-wise weak divergence (trial space ByDim)
  std::vector<mfem::DenseMatrix> weak_divergence_;
  /// element-wise inverse mass matrix (trial space ByDim)
  std::vector<mfem::DenseMatrix> inverse_mass_;
  /// global maximum characteristic speed. Updated by form integrators
  mutable mfem::real_t max_characteristic_speed_;
  /// Flux function \f$F(U)\f$ 
  std::shared_ptr<mfem::FluxFunction> flux_function_;
  /// Numerical flux function for computing \f$\hat{F} \cdot n\f$
  std::shared_ptr<mfem::NumericalFlux> numerical_flux_function_;
  /// Element integration form. Should contain ComputeFlux
  std::shared_ptr<mfem::HyperbolicFormIntegrator> hyperbolic_form_integrator_;
  /// Number of equations
  int num_equations_;
  /// ghost boundary condition setters
  std::vector<std::shared_ptr<DGGhostBC>> ghost_bcs_;
  /// ghost boundary condition integrators
  std::vector<std::shared_ptr<mfem::NonlinearFormIntegrator>> ghost_bc_integrators_;

  /// Compute element-wise inverse mass matrix
  void computeInvMass_();
  /// Compute element-wise weak-divergence matrix
  void computeWeakDivergence_();

};

} // namespace mfpic
