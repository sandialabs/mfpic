#include <libmfpic/ElectrostaticFieldOperationsWithBoltzmannElectrons.hpp>

namespace mfpic {

namespace {

/// Integrates the weak form of the Boltzmann electrons modifier term.
class BoltzmannElectronIntegrator : public mfem::NonlinearFormIntegrator {
public:
  /**
   * @brief Ctor.
   *
   * @param[in] reference_number_density Reference number density.
   * @param[in] temperature              Temperature.
   */
  BoltzmannElectronIntegrator(double reference_number_density, double temperature) :
    reference_number_density_(reference_number_density),
    temperature_(temperature)
  {}

  /**
   * @brief Assembles the weak form of the Boltzmann electron charge density.
   *
   * \f[
   *   e n_0 \int_\Omega \psi_i \exp \left( \frac{e \varphi_h}{k T} \right) \, \mathrm{d} \vec{x}
   * \f]
   *
   * @param[in]  element                                 Finite element in which to integrate.
   * @param[in]  element_transformation                  Physical-to-reference element transformation.
   * @param[in]  local_potential                         Potential vector on element-local dofs.
   * @param[out] local_boltzmann_electron_charge_density Electron charge density on element-local dofs.
   */
  virtual void AssembleElementVector(
    const mfem::FiniteElement& element,
    mfem::ElementTransformation& element_transformation,
    const mfem::Vector& local_potential,
    mfem::Vector& local_boltzmann_electron_charge_density
  ) override {
    const int integration_order = 2 * element.GetOrder() + element_transformation.OrderGrad(&element) + 1;
    IntRule = &mfem::IntRules.Get(element.GetGeomType(), integration_order);

    const int num_dofs = element.GetDof();
    local_boltzmann_electron_charge_density.SetSize(num_dofs);
    local_boltzmann_electron_charge_density = 0.0;

    mfem::Vector shape_functions_at_integration_point(num_dofs);
    for (int idof = 0; idof < IntRule->GetNPoints(); idof++) {
      const mfem::IntegrationPoint& integration_point = IntRule->IntPoint(idof);
      element_transformation.SetIntPoint(&integration_point);
      const double weight = element_transformation.Weight() * integration_point.weight;

      element.CalcPhysShape(element_transformation, shape_functions_at_integration_point);

      const double potential_at_integration_point = shape_functions_at_integration_point * local_potential;
      const double boltzmann_electron_charge_density_at_integration_point = reference_number_density_ * constants::elementary_charge * std::exp(
        constants::elementary_charge * potential_at_integration_point /
        (constants::boltzmann_constant * temperature_)
      );

      local_boltzmann_electron_charge_density.Add(
        weight * boltzmann_electron_charge_density_at_integration_point,
        shape_functions_at_integration_point
      );
    }
  }


  /**
   * @brief Assembles the weak form of the Boltzmann electron charge density Jacobian.
   *
   * \f[
   *   \frac{e^2 n_0}{kT} \int_\Omega \psi_i \psi_j \exp \left( \frac{e \varphi_h}{k T} \right) \, \mathrm{d} \vec{x}
   * \f]
   *
   * @param[in]  element                Finite element in which to integrate.
   * @param[in]  element_transformation Physical-to-reference element transformation.
   * @param[in]  local_potential        Potential vector on element-local dofs.
   * @param[out] local_jacobian         Jacobian of electron charge density on element-local dofs.
   */
  virtual void AssembleElementGrad(
    const mfem::FiniteElement& element,
    mfem::ElementTransformation& element_transformation,
    const mfem::Vector& local_potential,
    mfem::DenseMatrix& local_jacobian
  ) override {
    const int integration_order = 2 * element.GetOrder() + element_transformation.OrderGrad(&element) + 1;
    IntRule = &mfem::IntRules.Get(element.GetGeomType(), integration_order);

    const int num_dofs = element.GetDof();
    local_jacobian.SetSize(num_dofs, num_dofs);
    local_jacobian = 0.0;

    mfem::Vector shape_functions_at_integration_point(num_dofs);
    for (int idof = 0; idof < IntRule->GetNPoints(); idof++) {
      const mfem::IntegrationPoint& integration_point = IntRule->IntPoint(idof);
      element_transformation.SetIntPoint(&integration_point);
      const double weight = element_transformation.Weight() * integration_point.weight;

      element.CalcPhysShape(element_transformation, shape_functions_at_integration_point);

      const double potential_at_integration_point = shape_functions_at_integration_point * local_potential;
      const double jacobian_at_integration_point = reference_number_density_ * std::pow(constants::elementary_charge, 2.0)
        / (constants::boltzmann_constant * temperature_) * std::exp(
        constants::elementary_charge * potential_at_integration_point /
        (constants::boltzmann_constant * temperature_)
        );

      mfem::AddMult_a_VVt(weight * jacobian_at_integration_point, shape_functions_at_integration_point, local_jacobian);
    }
  }

private:
  /// Reference number density.
  const double reference_number_density_;

  /// Temperature.
  const double temperature_;
};

} // namespace

ElectrostaticFieldOperationsWithBoltzmannElectrons::ElectrostaticFieldOperationsWithBoltzmannElectrons(
  Discretization& electrostatic_discretization,
  std::unique_ptr<DirichletBoundaryConditions> dirichlet_boundary_conditions,
  double electron_reference_number_density,
  double electron_temperature,
  double nonlinear_solver_relative_tolerance,
  int nonlinear_solver_max_iterations
) :
  ElectrostaticFieldOperations(electrostatic_discretization, std::move(dirichlet_boundary_conditions)),
  electrostatic_nonlinear_form_(&electrostatic_discretization.getFeSpace())
{
  electrostatic_nonlinear_form_.AddDomainIntegrator(new mfem::DiffusionIntegrator(permittivity_));
  electrostatic_nonlinear_form_.AddDomainIntegrator(new BoltzmannElectronIntegrator(
    electron_reference_number_density,
    electron_temperature
  ));
  mfem::Array<int> dirichlet_dof_indices = dirichlet_boundary_conditions_->getDirichletBoundaryDofIndices();
  electrostatic_nonlinear_form_.SetEssentialTrueDofs(dirichlet_dof_indices);

  newton_solver_.SetOperator(electrostatic_nonlinear_form_);
  newton_solver_.SetSolver(cg_linear_solver_);
  newton_solver_.SetRelTol(nonlinear_solver_relative_tolerance);
  newton_solver_.SetMaxIter(nonlinear_solver_max_iterations);
}

void ElectrostaticFieldOperationsWithBoltzmannElectrons::fieldSolve(
  ElectrostaticFieldState& field_state,
  const IntegratedCharge& charge_state
) {
  mfem::Vector integrated_charge_vector = charge_state.getIntegratedCharge();
  const mfem::Array<int> dirichlet_dof_indices = dirichlet_boundary_conditions_->getDirichletBoundaryDofIndices();
  integrated_charge_vector.SetSubVector(dirichlet_dof_indices, 0.0);

  mfem::GridFunction& potential = field_state.getPotential();
  dirichlet_boundary_conditions_->applyBoundaryConditions(potential);

  newton_solver_.Mult(integrated_charge_vector, potential);
}

} // namespace mfpic
