#include <libmfpic/Constants.hpp>
#include <libmfpic/ElectrostaticFieldOperations.hpp>
#include <libmfpic/MFEMHelpers.hpp>

#include <mfem/fem/bilininteg.hpp>
#include <mfem/fem/gridfunc.hpp>

namespace mfpic {

ElectrostaticFieldOperations::ElectrostaticFieldOperations(
  Discretization& electrostatic_discretization,
  std::unique_ptr<DirichletBoundaryConditions> dirichlet_boundary_conditions)
  : dirichlet_boundary_conditions_(std::move(dirichlet_boundary_conditions))
  , electrostatic_bilinear_form_(&electrostatic_discretization.getFeSpace())
  , mass_form_(&electrostatic_discretization.getFeSpace())
  , scalar_dg_discretization_(electrostatic_discretization.getFeSpace().GetMesh(), 0, FETypes::DG)
  , weak_derivative_bilinear_form_(&scalar_dg_discretization_.getFeSpace(), &electrostatic_discretization.getFeSpace())
{
  mfem::ConstantCoefficient permittivity(constants::permittivity);
  electrostatic_bilinear_form_.AddDomainIntegrator(new mfem::DiffusionIntegrator(permittivity));
  electrostatic_bilinear_form_.Assemble();

  mass_form_.AddDomainIntegrator(new mfem::LumpedIntegrator(new mfem::MassIntegrator));
  mass_form_.Assemble();

  weak_derivative_bilinear_form_.AddDomainIntegrator(new mfem::MixedScalarWeakDerivativeIntegrator(permittivity));
  weak_derivative_bilinear_form_.Assemble();
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
  double field_energy;
  if (potential.Norml2() > 0) {
    field_energy = 0.5 * electrostatic_bilinear_form_.InnerProduct(potential, potential);
  } else { // if potential and e_field are both zero then this computation will give zero field energy
    const mfem::GridFunction& e_field = field_state.getEFieldGridFunction();

    auto e_field_coefficient = std::make_unique<mfem::VectorGridFunctionCoefficient>(&e_field);
    auto field_energy_function = [](const mfem::Vector& E){ return 0.5 * constants::permittivity * (E * E); };
    TransformedVectorCoefficient field_energy_coefficient(std::move(e_field_coefficient), field_energy_function);

    const mfem::FiniteElementSpace* e_field_fe_space = e_field.FESpace();
    mfem::Mesh& mesh = *(e_field_fe_space->GetMesh());
    Discretization identity_test_function_discretization = getIdentityTestFunctionDiscretization(mesh);
    mfem::LinearForm field_energy_by_cell(&identity_test_function_discretization.getFeSpace());
    field_energy_by_cell.AddDomainIntegrator(new mfem::DomainLFIntegrator(field_energy_coefficient));
    field_energy_by_cell.Assemble();

    field_energy = field_energy_by_cell.Sum();
  }
  return field_energy;
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

mfem::Vector ElectrostaticFieldOperations::computeIntegratedGhostCharge2(
  const ElectrostaticFieldState& es_field_state,
  const IntegratedCharge& integrated_charge)
{
  const mfem::GridFunction& e_field = es_field_state.getEFieldGridFunction();
  const mfem::Vector integrated_charge_vector = integrated_charge.getIntegratedCharge();

  mfem::Vector e_field_copy(e_field);
  const mfem::GridFunction e_field_x(&scalar_dg_discretization_.getFeSpace(), e_field_copy);

  mfem::Vector integrated_ghost_charge(integrated_charge_vector.Size());
  weak_derivative_bilinear_form_.Mult(e_field_x, integrated_ghost_charge);
  integrated_ghost_charge.Add(-1, integrated_charge_vector);

  return integrated_ghost_charge;
}

mfem::GridFunction ElectrostaticFieldOperations::computeGhostChargeDensity2(
  const ElectrostaticFieldState& es_field_state,
  const IntegratedCharge& integrated_charge)
{
  mfem::Vector integrated_ghost_charge = computeIntegratedGhostCharge2(es_field_state, integrated_charge);

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

mfem::GridFunction ElectrostaticFieldOperations::computeWeakDivE(
  const ElectrostaticFieldState& es_field_state)
{
  const mfem::GridFunction& e_field = es_field_state.getEFieldGridFunction();

  mfem::Vector e_field_copy(e_field);
  const mfem::GridFunction e_field_x(&scalar_dg_discretization_.getFeSpace(), e_field_copy);

  mfem::Vector weak_div_e_integrated(e_field_x.Size());
  weak_derivative_bilinear_form_.Mult(e_field_x, weak_div_e_integrated);

  mfem::GridFunction weak_div_e(mass_form_.FESpace());

  mfem::Vector solution_vector;
  mfem::Vector rhs_vector;
  mfem::Array<int> no_bcs;
  mass_form_.FormLinearSystem(
    no_bcs,
    weak_div_e,
    weak_div_e_integrated,
    mass_matrix_,
    solution_vector,
    rhs_vector);

  cg_linear_solver_.solve(mass_matrix_, solution_vector, rhs_vector);
  mass_form_.RecoverFEMSolution(solution_vector, weak_div_e_integrated, weak_div_e);

  return weak_div_e;

}

ElectrostaticFieldOperations::~ElectrostaticFieldOperations() = default;

}
