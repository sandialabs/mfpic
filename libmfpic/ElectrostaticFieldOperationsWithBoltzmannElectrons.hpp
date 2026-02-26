#pragma once

#include <libmfpic/ElectrostaticFieldOperations.hpp>

namespace mfpic {

class BoltzmannElectronIntegrator : public mfem::NonlinearFormIntegrator {
public:
  BoltzmannElectronIntegrator(double reference_number_density, double temperature);

  virtual void AssembleElementVector(
    const mfem::FiniteElement& element,
    mfem::ElementTransformation& element_transformation,
    const mfem::Vector& local_potential,
    mfem::Vector& local_boltzmann_electron_charge_density
  );

  virtual void AssembleElementGrad(
    const mfem::FiniteElement& element,
    mfem::ElementTransformation& element_transformation,
    const mfem::Vector& local_potential,
    mfem::DenseMatrix& local_jacobian
  );

private:
  const double reference_number_density_;

  const double temperature_;
};

class ElectrostaticFieldOperationsWithBoltzmannElectrons : public ElectrostaticFieldOperations {
public:
  ElectrostaticFieldOperationsWithBoltzmannElectrons(
    Discretization& electrostatic_discretization,
    std::unique_ptr<DirichletBoundaryConditions> dirichlet_boundary_conditions
  );

  virtual void fieldSolve(ElectrostaticFieldState& field_state, const IntegratedCharge& charge_state) override;

private:
  mfem::ConstantCoefficient permittivity_ = mfem::ConstantCoefficient(constants::permittivity);

  mfem::NonlinearForm electrostatic_nonlinear_form_;
};

} // namespace mfpic
