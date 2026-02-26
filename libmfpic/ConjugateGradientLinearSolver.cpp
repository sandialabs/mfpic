#include <libmfpic/ConjugateGradientLinearSolver.hpp>
#include <libmfpic/Errors.hpp>

namespace mfpic {

ConjugateGradientLinearSolver::ConjugateGradientLinearSolver(
  const double relative_tolerance,
  const double absolute_tolerance,
  const int maximum_number_of_iterations)
  : solver_log_("solver.log")
{
  SetPreconditioner(gauss_seidel_preconditioner_);
  SetRelTol(relative_tolerance);
  SetAbsTol(absolute_tolerance);
  SetMaxIter(maximum_number_of_iterations);
  SetPrintLevel(mfem::IterativeSolver::PrintLevel().Warnings());
}

void ConjugateGradientLinearSolver::solve(
  const mfem::SparseMatrix& matrix_to_solve,
  mfem::Vector& sol_vector,
  const mfem::Vector& rhs_vector)
{
  SetOperator(matrix_to_solve);
  Mult(rhs_vector, sol_vector);

  const int num_iterations = GetNumIterations();
  const double initial_residual_norm = GetInitialNorm();
  const double final_residual_norm = GetFinalNorm();
  const double final_residual_relative_norm = GetFinalRelNorm();

  std::ostringstream summary;
  summary << "Linear Solve Summary:\n";
  summary << "  Number of Iterations = " << num_iterations << "\n";
  summary << "  Initial Residual Norm = " << initial_residual_norm << "\n";
  summary << "  Final Residual Relative Norm = " << final_residual_relative_norm << "\n";
  summary << "  Final Residual Norm = " << final_residual_norm << "\n";
  summary << "  Relative Tolerance = " << rel_tol << "\n";
  summary << "  Absolute Tolerance = " << abs_tol << "\n";
  solver_log_ << summary.str();

  const bool solve_converged = GetConverged();
  if (not solve_converged) {
    std::ostringstream error_message;
    error_message << "Conjugate Gradient Linear Solver failed to converge.\n";
    error_message << summary.str();

    errorWithUserMessage(error_message.str());
  }
}

}