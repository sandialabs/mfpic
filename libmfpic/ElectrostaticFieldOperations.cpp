#include <libmfpic/Constants.hpp>
#include <libmfpic/ElectrostaticFieldOperations.hpp>
#include <mfem/fem/bilininteg.hpp>
#include <mfem/fem/gridfunc.hpp>

namespace mfpic {

ElectrostaticFieldOperations::ElectrostaticFieldOperations(
  Discretization& electrostatic_discretization,
  std::unique_ptr<DirichletBoundaryConditions> dirichlet_boundary_conditions)
  : dirichlet_boundary_conditions_(std::move(dirichlet_boundary_conditions))
  , electrostatic_bilinear_form_(&electrostatic_discretization.getFeSpace())
  , mass_form_(&electrostatic_discretization.getFeSpace())
{
  mfem::ConstantCoefficient permittivity(constants::permittivity);
  electrostatic_bilinear_form_.AddDomainIntegrator(new mfem::DiffusionIntegrator(permittivity));
  electrostatic_bilinear_form_.Assemble();

  mass_form_.AddDomainIntegrator(new mfem::MassIntegrator);
  mass_form_.Assemble();
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

  mfem::GridFunction charge_from_integrated_charge_grid_function(potential_copy.FESpace());
  mfem::GridFunction charge_from_potential_grid_function(potential_copy.FESpace());

  auto compute_charge = [&](mfem::Vector& rhs, mfem::GridFunction& result) {
    mfem::Vector x;
    mfem::Vector b;
    mfem::Array<int> no_bcs;
    mass_form_.FormLinearSystem(no_bcs, 
                                result, 
                                rhs, 
                                mass_matrix_, 
                                x, 
                                b);
    cg_linear_solver_.solve(mass_matrix_, x, b);
    mass_form_.RecoverFEMSolution(x, rhs, result);
  };

  auto copy_boundary_dofs = [&](const mfem::Vector& source, mfem::Vector& target) {
    mfem::Array<int> boundary_attribute_flags;
    boundary_attribute_flags.SetSize(potential_copy.FESpace()->GetMesh()->bdr_attributes.Max());
    boundary_attribute_flags = true;
    mfem::Array<int> all_boundary_dof_indices;
    potential_copy.FESpace()->GetEssentialTrueDofs(boundary_attribute_flags, all_boundary_dof_indices);

    for (int i = 0; i < all_boundary_dof_indices.Size(); ++i)
      target(all_boundary_dof_indices[i]) = source(all_boundary_dof_indices[i]);
  };

  mfem::Vector charge_vector(potential_copy.Size());
  electrostatic_bilinear_form_.Mult(potential_copy, charge_vector);

  mfem::Vector integrated_charge_vector = integrated_charge.getIntegratedCharge();
  copy_boundary_dofs(integrated_charge_vector, charge_vector); // TODO BWR ADD NOTE ABOUT THIS IF WE KEEP THIS

  compute_charge(charge_vector, charge_from_potential_grid_function);
  compute_charge(integrated_charge_vector, charge_from_integrated_charge_grid_function);

  mfem::GridFunction diff = charge_from_integrated_charge_grid_function;
  diff -= charge_from_potential_grid_function;

  return diff;

}

ElectrostaticFieldOperations::~ElectrostaticFieldOperations() = default;

}
