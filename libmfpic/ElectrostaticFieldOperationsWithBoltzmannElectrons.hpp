#pragma once

#include "libmfpic/Constants.hpp"
#include <libmfpic/ElectrostaticFieldOperations.hpp>

namespace mfpic {

class BoltzmannElectronIntegrator : public mfem::NonlinearFormIntegrator {
public:
  virtual void AssembleElementVector(
    const mfem::FiniteElement& element,
    mfem::ElementTransformation& element_transformation,
    const mfem::Vector& local_potential,
    mfem::Vector& local_boltzmann_electron_charge_density
  ) {
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

  virtual void AssembleElementGrad(
    const mfem::FiniteElement& element,
    mfem::ElementTransformation& element_transformation,
    const mfem::Vector& local_potential,
    mfem::DenseMatrix& local_jacobian
  ) {
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
  const double reference_number_density_ = 3.54803829431936e+18;
  // const double reference_number_density_ = 7.05e18;

  const double temperature_ = 10.0 * constants::elementary_charge / constants::boltzmann_constant;
};

class ElectrostaticFieldOperationsWithBoltzmannElectrons : public ElectrostaticFieldOperations {
public:
  ElectrostaticFieldOperationsWithBoltzmannElectrons(
    Discretization& electrostatic_discretization,
    std::unique_ptr<DirichletBoundaryConditions> dirichlet_boundary_conditions
  ) :
    ElectrostaticFieldOperations(electrostatic_discretization, std::move(dirichlet_boundary_conditions)),
    electrostatic_nonlinear_form_(&electrostatic_discretization.getFeSpace())
  {
    electrostatic_nonlinear_form_.AddDomainIntegrator(new mfem::DiffusionIntegrator(permittivity_));
    electrostatic_nonlinear_form_.AddDomainIntegrator(new BoltzmannElectronIntegrator());
    mfem::Array<int> dirichlet_dof_indices = dirichlet_boundary_conditions_->getDirichletBoundaryDofIndices();
    electrostatic_nonlinear_form_.SetEssentialTrueDofs(dirichlet_dof_indices);
  }

  virtual void fieldSolve(ElectrostaticFieldState& field_state, const IntegratedCharge& charge_state) override {
    const mfem::Vector integrated_charge_vector = charge_state.getIntegratedCharge();
    mfem::GridFunction& potential = field_state.getPotential();

    mfem::CGSolver linear_solver;
    mfem::NewtonSolver newton_solver;
    newton_solver.SetOperator(electrostatic_nonlinear_form_);
    newton_solver.SetSolver(linear_solver);
    newton_solver.SetRelTol(1.0e-8);
    newton_solver.SetMaxIter(100);
    newton_solver.Mult(integrated_charge_vector, potential);
  }

private:
  mfem::ConstantCoefficient permittivity_ = mfem::ConstantCoefficient(constants::permittivity);

  mfem::NonlinearForm electrostatic_nonlinear_form_;
};

} // namespace mfpic
