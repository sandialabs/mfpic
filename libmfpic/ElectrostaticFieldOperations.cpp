#include <libmfpic/Constants.hpp>
#include <libmfpic/ElectrostaticFieldOperations.hpp>

#include <mfem.hpp>

namespace mfpic {

ElectrostaticFieldOperations::ElectrostaticFieldOperations(
  Discretization& electrostatic_discretization,
  std::unique_ptr<DirichletBoundaryConditions> dirichlet_boundary_conditions)
  : dirichlet_boundary_conditions_(std::move(dirichlet_boundary_conditions))
  , electrostatic_bilinear_form_(&electrostatic_discretization.getFeSpace())
  , mass_form_(&electrostatic_discretization.getFeSpace())
  , null_space_(&electrostatic_discretization.getFeSpace())
{
  mfem::ConstantCoefficient permittivity(constants::permittivity);
  electrostatic_bilinear_form_.AddDomainIntegrator(new mfem::DiffusionIntegrator(permittivity));
  electrostatic_bilinear_form_.Assemble();

  mass_form_.AddDomainIntegrator(new mfem::LumpedIntegrator(new mfem::MassIntegrator));
  mass_form_.Assemble();

  auto one = mfem::ConstantCoefficient(1.0);
  null_space_.ProjectCoefficient(one);
}

void ElectrostaticFieldOperations::fieldSolve(ElectrostaticFieldState& field_state, const IntegratedCharge& charge_state) {
  mfem::GridFunction potential = field_state.getPotential();

  dirichlet_boundary_conditions_->applyBoundaryConditions(potential);

  mfem::Vector integrated_charge_vector = charge_state.getIntegratedCharge();

  mfem::Array<int> dirichlet_dof_indices = dirichlet_boundary_conditions_->getDirichletBoundaryDofIndices();
  if (dirichlet_dof_indices.IsEmpty()) {
    enforceCompatibilityOnIntegratedCharge(integrated_charge_vector);
  }
  
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

mfem::Vector ElectrostaticFieldOperations::computeIntegratedGhostCharge(
  const ElectrostaticFieldState& es_field_state,
  const IntegratedCharge& integrated_charge) const
{
  const mfem::GridFunction& potential = es_field_state.getPotential();
  const mfem::Vector integrated_charge_vector = integrated_charge.getIntegratedCharge();

  mfem::Vector integrated_ghost_charge(potential.Size());

  electrostatic_bilinear_form_.FullMult(potential, integrated_ghost_charge);
  integrated_ghost_charge.Add(-1., integrated_charge_vector);

  return integrated_ghost_charge;
}

mfem::GridFunction ElectrostaticFieldOperations::computeGhostChargeDensity(
  const ElectrostaticFieldState& es_field_state,
  const IntegratedCharge& integrated_charge)
{
  mfem::Vector integrated_ghost_charge = computeIntegratedGhostCharge(es_field_state, integrated_charge);

  mfem::GridFunction ghost_charge_density(mass_form_.FESpace());

  mfem::Vector solution_vector;
  mfem::Vector rhs_vector;
  mfem::Array<int> no_bcs;
  mass_form_.FormLinearSystem(
    no_bcs,
    ghost_charge_density,
    integrated_ghost_charge,
    mass_matrix_,
    solution_vector,
    rhs_vector);

  cg_linear_solver_.solve(mass_matrix_, solution_vector, rhs_vector);
  mass_form_.RecoverFEMSolution(solution_vector, integrated_ghost_charge, ghost_charge_density);

  return ghost_charge_density;
}

void ElectrostaticFieldOperations::enforceCompatibilityOnIntegratedCharge(mfem::Vector& integrated_charge_vector) const 
{
  const double rhs_dot_null = null_space_ * integrated_charge_vector / (null_space_ * null_space_);
  integrated_charge_vector.Add(-rhs_dot_null, null_space_);
}

ElectrostaticFieldOperations::~ElectrostaticFieldOperations() = default;

}
