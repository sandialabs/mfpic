#pragma once

#include <mfem/mfem.hpp>

namespace mfpic {

class ConjugateGradientLinearSolver {
public:
  /**
   * @brief Construct a new Conjugate Gradient Linear Solver object
   * 
   * Always uses jacobi preconditioner
   * 
   * @param relative_tolerance - iterative solve stops if residual_norm / initial_residual_norm < relative_tolerance
   * @param absolute_tolerance - iterative solve stops if residual_norm < absolute_tolerance
   * @param maximum_number_of_iterations - maximum number of iterations for iterative solver
   */
  ConjugateGradientLinearSolver(
    const double relative_tolerance = 1e-12,
    const double absolute_tolerance = 0,
    const int maximum_number_of_iterations = 1000);

  /**
   * @brief Solve the linear system, throws error if solver fails to converge
   * 
   * @note It is possible that linear solver will fail but still report convergence. In that case an error message will be output 
   * to screen but code will continue to run.
   *
   * @param matrix_to_solve - matrix to invert
   * @param solution_vector - vector to store solution of linear system, provides initial guess for iterative solver
   * @param rhs_vector - vector that stores right hand side of linear system
   */
  void solve(const mfem::SparseMatrix& matrix_to_solve, mfem::Vector& solution_vector, const mfem::Vector& rhs_vector);

private:
  mfem::GSSmoother gauss_seidel_preconditioner_;
  mfem::CGSolver mfem_cg_solver_;

  double relative_tolerance_;
  double absolute_tolerance_;
  int maximum_number_of_iterations_;

  std::ofstream solver_log_;
};

}