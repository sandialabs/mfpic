#pragma once

#include <libmfpic/Constants.hpp>
#include <libmfpic/ElectrostaticFieldOperations.hpp>

namespace mfpic {

/// Electrostatic field solver incorporating Boltzmann electrons.
class ElectrostaticFieldOperationsWithBoltzmannElectrons : public ElectrostaticFieldOperations {
public:
  /**
   * @brief Ctor.
   *
   * @param[in] electrostatic_discretization        Potential finite element space to determine bilinear form
   *                                                This object is passed in non const because mfem::BilinearForm requires a non const FiniteElementSpace, however the
   *                                                BilinearForm shouldn't change any meaningful data.
   * @param[in] dirichlet_boundary_conditions       Dirichlet boundary conditions for field solve, transfers ownership to this object.
   * @param[in] electron_reference_number_density   Reference number density of the electrons.
   * @param[in] electron_temperature                Temperature of the electrons.
   * @param[in] nonlinear_solver_relative_tolerance Relative tolerance for the nonlinear solver.
   * @param[in] nonlinear_solver_max_iterations     Maximum iterations for the nonlinear solver.
   */
  ElectrostaticFieldOperationsWithBoltzmannElectrons(
    Discretization& electrostatic_discretization,
    std::unique_ptr<DirichletBoundaryConditions> dirichlet_boundary_conditions,
    double electron_reference_number_density,
    double electron_temperature,
    double nonlinear_solver_relative_tolerance = 1.0e-8,
    int nonlinear_solver_max_iterations = 100
  );

  /**
   * @brief solve for the field state that satisfies the electrostatic model with a given charge state
   *
   * @note this method is non const because mfem::BilinearForm must be non const to FormLinearSystem.
   *  However no meaningful internal data should be changed.
   *
   * @param field_state - where to store new field state, also initial guess for linear solver
   * @param charge_state - integrated charge for right hand side of system, \int_{\Omega} \rho \phi dV
   */
  virtual void fieldSolve(ElectrostaticFieldState& field_state, const IntegratedCharge& charge_state) override;

private:
  /// Vacuum permittivity.
  mfem::ConstantCoefficient permittivity_ = mfem::ConstantCoefficient(constants::permittivity);

  /// Nonlinear form consisting of the discrete Laplacian and Boltzmann electron modifier.
  mfem::NonlinearForm electrostatic_nonlinear_form_;

  /// Newton solver used for the nonlinear field solve.
  mfem::NewtonSolver newton_solver_;
};

} // namespace mfpic
