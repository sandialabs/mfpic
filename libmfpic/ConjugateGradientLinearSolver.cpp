#include <libmfpic/ConjugateGradientLinearSolver.hpp>
#include <libmfpic/Errors.hpp>

namespace mfpic {

ConjugateGradientLinearSolver::ConjugateGradientLinearSolver(
  const double relative_tolerance,
  const double absolute_tolerance,
  const int maximum_number_of_iterations)
  : relative_tolerance_(relative_tolerance)
  , absolute_tolerance_(absolute_tolerance)
  , maximum_number_of_iterations_(maximum_number_of_iterations)
  , solver_log_("solver.log")
{
  mfem_cg_solver_.SetPreconditioner(gauss_seidel_preconditioner_);
  mfem_cg_solver_.SetRelTol(relative_tolerance_);
  mfem_cg_solver_.SetAbsTol(absolute_tolerance_);
  mfem_cg_solver_.SetMaxIter(maximum_number_of_iterations_);
  mfem_cg_solver_.SetPrintLevel(mfem::IterativeSolver::PrintLevel().Warnings());
}

void ConjugateGradientLinearSolver::solve(
  const mfem::SparseMatrix& matrix_to_solve,
  mfem::Vector& sol_vector,
  const mfem::Vector& rhs_vector)
{
  mfem_cg_solver_.SetOperator(matrix_to_solve);
  mfem_cg_solver_.Mult(rhs_vector, sol_vector);

  const int num_iterations = mfem_cg_solver_.GetNumIterations();
  const double initial_residual_norm = mfem_cg_solver_.GetInitialNorm();
  const double final_residual_norm = mfem_cg_solver_.GetFinalNorm();
  const double final_residual_relative_norm = mfem_cg_solver_.GetFinalRelNorm();

  std::ostringstream summary;
  summary << "Linear Solve Summary:\n";
  summary << "  Number of Iterations = " << num_iterations << "\n";
  summary << "  Initial Residual Norm = " << initial_residual_norm << "\n";
  summary << "  Final Residual Relative Norm = " << final_residual_relative_norm << "\n";
  summary << "  Final Residual Norm = " << final_residual_norm << "\n";
  summary << "  Relative Tolerance = " << relative_tolerance_ << "\n";
  summary << "  Absolute Tolerance = " << absolute_tolerance_ << "\n";
  solver_log_ << summary.str();

  const bool solve_converged = mfem_cg_solver_.GetConverged();
  if (not solve_converged) {
    std::ostringstream error_message;
    error_message << "Conjugate Gradient Linear Solver failed to converge.\n";
    error_message << summary.str();

    errorWithUserMessage(error_message.str());
  }
}

}