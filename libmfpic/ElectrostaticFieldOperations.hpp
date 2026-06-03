#pragma once

#include <libmfpic/ConjugateGradientLinearSolver.hpp>
#include <libmfpic/DirichletBoundaryConditions.hpp>
#include <libmfpic/IntegratedCharge.hpp>
#include <libmfpic/ElectrostaticFieldState.hpp>

#include <mfem/fem/bilinearform.hpp>
#include <mfem/mfem.hpp>

namespace mfpic {

class ElectrostaticFieldOperations {
public:
  /**
   * @brief Construct a new Electrostatic Field Operations object
   * 
   * @param electrostatic_discretization - potential finite element space to determine bilinear form
   *  This object is passed in non const because mfem::BilinearForm requires a non const FiniteElementSpace, however the
   *  BilinearForm shouldn't change any meaningful data.
   * @param dirichlet_boundary_conditions - dirichlet boundary conditions for field solve, transfers ownership to this object
   */
  ElectrostaticFieldOperations(
    Discretization& electrostatic_discretization,
    std::unique_ptr<DirichletBoundaryConditions> dirichlet_boundary_conditions);

  /**
   * @brief solve for the field state that satisfies the electrostatic model with a given charge state
   * 
   * @note this method is non const because mfem::BilinearForm must be non const to FormLinearSystem.
   *  However no meaningful internal data should be changed.
   * 
   * @param field_state - where to store new field state, also initial guess for linear solver
   * @param charge_state - integrated charge for right hand side of system, \int_{\Omega} \rho \psi dV
   */
  virtual void fieldSolve(ElectrostaticFieldState& field_state, const IntegratedCharge& charge_state);

  /**
   * @brief Compute the energy of the given field state.
   *
   * @param[in] field_state Field for which to compute the energy.
   *
   * @return The electrostatic field energy in the domain.
   */
  double fieldEnergy(const ElectrostaticFieldState& field_state) const;

  /**
   * @brief Compute the integrated ghost charge \f$-L_{\eps} \phi - \rho_I\f$
   * where \f$L_{\eps}\f$ is the Laplacian matrix weighted by the permittivy and \f$\rho_I\f$ is the
   * given integrated charge
   *
   * @param field_state \ref ElectrostaticFieldState supplying the potential \f$\phi\f$
   * @param integrated_charge \ref IntegratedCharge \f$\rho_I\f$
   *
   * @returns Integrated ghost charge
   */

  mfem::Vector computeIntegratedGhostCharge(
    const ElectrostaticFieldState& field_state,
    const IntegratedCharge& integrated_charge) const;

  /**
   * @brief Compute the ghost charge density \f$M^{-1}\left(-L_{\eps} \phi - \rho_I\right)\f$
   * where \f$L_{\eps}\f$ is the Laplacian matrix weighted by the permittivy, \f$\rho_I\f$ is the
   * given integrated charge, and \f$M\f$ is the mass matrix
   *
   * @param field_state \ref ElectrostaticFieldState supplying the potential \f$\phi\f$
   * @param integrated_charge \ref IntegratedCharge \f$\rho_I\f$
   *
   * @returns Ghost charge density
   */
  mfem::GridFunction computeGhostChargeDensity(
    const ElectrostaticFieldState& field_state,
    const IntegratedCharge& integrated_charge);

  /**
   * @brief Remove the component of the integrated charge in the direction of the null space of the Laplacian.
   *
   * @param[in,out] integrated_charge_vector mfem::Vector of the integrated charge
   */
  void enforceCompatabilityOnIntegratedCharge(mfem::Vector& integrated_charge_vector) const;

  /// Dtor.
  virtual ~ElectrostaticFieldOperations();

protected:
  /// Dirichlet boundary conditions.
  std::unique_ptr<DirichletBoundaryConditions> dirichlet_boundary_conditions_;

  /// conjugate gradient solver to invert negative_eps_laplace_matrix_
  ConjugateGradientLinearSolver cg_linear_solver_;

private:
  /// bilinear form for electrostatic system, -eps Laplace or eps (grad potential, grad phi)
  mfem::BilinearForm electrostatic_bilinear_form_;

  /// -eps Laplace operator assembled into matrix
  mfem::SparseMatrix negative_eps_laplace_matrix_;

  /// bilinear mass form for computing ghost charge
  mfem::BilinearForm mass_form_;

  /// mass operator assembled into matrix
  mfem::SparseMatrix mass_matrix_;

  /// grid function representing the null space of the Laplacian
  mfem::GridFunction null_space_;
};

} // namespace mfpic
