#include <libmfpic/Constants.hpp>
#include <libmfpic/ElectrostaticFieldOperations.hpp>

namespace mfpic {

ElectrostaticFieldOperations::ElectrostaticFieldOperations(
  Discretization& electrostatic_discretization,
  std::unique_ptr<DirichletBoundaryConditions> dirichlet_boundary_conditions)
  : dirichlet_boundary_conditions_(std::move(dirichlet_boundary_conditions))
  , electrostatic_bilinear_form_(&electrostatic_discretization.getFeSpace())
{
  mfem::ConstantCoefficient permittivity(constants::permittivity);
  electrostatic_bilinear_form_.AddDomainIntegrator(new mfem::DiffusionIntegrator(permittivity));
  electrostatic_bilinear_form_.Assemble();
}

void ElectrostaticFieldOperations::fieldSolve(ElectrostaticFieldState& field_state, const IntegratedCharge& charge_state) {
  mfem::GridFunction potential = field_state.getPotential();

  dirichlet_boundary_conditions_->applyBoundaryConditions(potential);

  mfem::Vector integrated_charge_vector = charge_state.getIntegratedCharge();

  mfem::Array<int> dirichlet_dof_indices = dirichlet_boundary_conditions_->getDirichletBoundaryDofIndices();
  mfem::Vector solution_vector;
  mfem::Vector rhs_vector;
  electrostatic_bilinear_form_.FormLinearSystem(
    dirichlet_dof_indices,
    potential,
    integrated_charge_vector,
    negative_eps_laplace_matrix_,
    solution_vector,
    rhs_vector);

  cg_linear_solver_.solve(negative_eps_laplace_matrix_, solution_vector, rhs_vector);

  electrostatic_bilinear_form_.RecoverFEMSolution(solution_vector, integrated_charge_vector, potential);

  field_state.setPotential(potential);
}

double ElectrostaticFieldOperations::fieldEnergy(const ElectrostaticFieldState& field_state) const {
  const mfem::GridFunction& potential = field_state.getPotential();
  return 0.5 * electrostatic_bilinear_form_.InnerProduct(potential, potential);
}

mfem::GridFunction ElectrostaticFieldOperations::chargeError(const ElectrostaticFieldState& field_state, const IntegratedCharge& integrated_charge) {
  mfem::GridFunction potential_copy = field_state.getPotential();

  mfem::GridFunction error(potential_copy.FESpace());

  mfem::Vector charge_deficit(potential_copy.Size());
  electrostatic_bilinear_form_.Mult(potential_copy, charge_deficit);
  charge_deficit -= integrated_charge.getIntegratedCharge();

  mfem::BilinearForm mass_form(potential_copy.FESpace());
  mass_form.AddDomainIntegrator(new mfem::MassIntegrator);
  mass_form.Assemble();
  mass_form.Finalize();

  mfem::Array<int> dirichlet_dof_indices; // TODO BWR not sure about BCs
  mfem::Vector solution_vector;
  mfem::Vector rhs_vector;
  mfem::SparseMatrix mass_matrix;
  mass_form.FormLinearSystem(dirichlet_dof_indices, 
                             error, 
                             charge_deficit, 
                             mass_matrix, 
                             solution_vector, 
                             rhs_vector);
  cg_linear_solver_.solve(mass_matrix, solution_vector, rhs_vector);
  mass_form.RecoverFEMSolution(solution_vector, charge_deficit, error);

  return error;

}

ElectrostaticFieldOperations::~ElectrostaticFieldOperations() = default;

}
