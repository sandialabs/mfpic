#pragma once

#include <libmfpic/Constants.hpp>
#include <libmfpic/ElectrostaticFieldOperations.hpp>

namespace mfpic {

class ElectrostaticFieldOperationsWithBoltzmannElectrons : public ElectrostaticFieldOperations {
public:
  ElectrostaticFieldOperationsWithBoltzmannElectrons(
    Discretization& electrostatic_discretization,
    std::unique_ptr<DirichletBoundaryConditions> dirichlet_boundary_conditions,
    double electron_reference_number_density,
    double electron_temperature,
    double nonlinear_solver_relative_tolerance = 1.0e-8,
    double nonlinear_solver_max_iterations = 100
  );

  virtual void fieldSolve(ElectrostaticFieldState& field_state, const IntegratedCharge& charge_state) override;

private:
  mfem::ConstantCoefficient permittivity_ = mfem::ConstantCoefficient(constants::permittivity);

  mfem::NonlinearForm electrostatic_nonlinear_form_;

  mfem::NewtonSolver newton_solver_;
};

} // namespace mfpic
