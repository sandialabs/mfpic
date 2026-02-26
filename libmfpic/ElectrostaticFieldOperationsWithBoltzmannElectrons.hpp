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
  );

  virtual void AssembleElementGrad(
    const mfem::FiniteElement& element,
    mfem::ElementTransformation& element_transformation,
    const mfem::Vector& local_potential,
    mfem::DenseMatrix& local_jacobian
  );

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
  );

  virtual void fieldSolve(ElectrostaticFieldState& field_state, const IntegratedCharge& charge_state) override;

private:
  mfem::ConstantCoefficient permittivity_ = mfem::ConstantCoefficient(constants::permittivity);

  mfem::NonlinearForm electrostatic_nonlinear_form_;
};

} // namespace mfpic
