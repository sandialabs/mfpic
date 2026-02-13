#include <libmfpic/DGAssembly.hpp>
#include <libmfpic/ElectromagneticFieldsEvaluator.hpp>
#include <libmfpic/Errors.hpp>
#include <mfem/fem/fespace.hpp>
#include <mfem/fem/gridfunc.hpp>
#include <mfem/fem/hyperbolic.hpp>
#include <mfem/linalg/densemat.hpp>

namespace mfpic {

  DGAssembly::DGAssembly(
    mfem::FiniteElementSpace &finite_element_space,
    std::unique_ptr<mfem::FluxFunction> &&flux_function,
    std::unique_ptr<mfem::NumericalFlux> &&numerical_flux_function) :
      dim_(finite_element_space.GetMesh()->SpaceDimension()),
      finite_element_space_(finite_element_space),
      flux_function_(std::move(flux_function)),
      numerical_flux_function_(std::move(numerical_flux_function)),
      hyperbolic_form_integrator_(std::make_unique<mfem::HyperbolicFormIntegrator>(*numerical_flux_function_)),
      num_equations_(hyperbolic_form_integrator_->num_equations)
    {
      if (finite_element_space_.GetOrdering() != mfem::Ordering::byNODES) {
        std::string error_message = "DGEulerOperator requires the finite element space be ordered 'byNODES'";
        errorWithDeveloperMessage(error_message);
      }
      this->computeInvMass_();
#ifndef MFEM_USE_MPI
      nonlinear_form_.reset(new mfem::NonlinearForm(&finite_element_space));
#else
    mfem::ParFiniteElementSpace *parallel_finite_element_space = dynamic_cast<mfem::ParFiniteElementSpace *>(&finite_element_space);
    if (parallel_finite_element_space)
    {
      nonlinear_form_.reset(new mfem::ParNonlinearForm(parallel_finite_element_space));
    }
    else
    {
      nonlinear_form_.reset(new mfem::NonlinearForm(&finite_element_space));
    }
#endif
    this->computeWeakDivergence_();
    nonlinear_form_->AddInteriorFaceIntegrator(hyperbolic_form_integrator_.get());
    nonlinear_form_->UseExternalIntegrators();
  }

  void DGAssembly::computeInvMass_()
  {
    mfem::InverseIntegrator inverse_mass_integrator(new mfem::MassIntegrator());

    inverse_mass_.resize(finite_element_space_.GetNE());
    for (int element=0; element<finite_element_space_.GetNE(); element++)
    {
      const int num_dof = finite_element_space_.GetFE(element)->GetDof();
      inverse_mass_[element].SetSize(num_dof);
      inverse_mass_integrator.AssembleElementMatrix(*finite_element_space_.GetFE(element),
                                                    *finite_element_space_.GetElementTransformation(element),
                                                    inverse_mass_[element]);
    }
  }

  void DGAssembly::computeWeakDivergence_()
  {
    mfem::TransposeIntegrator weak_divergence_integrator(new mfem::GradientIntegrator());
    mfem::DenseMatrix weak_divergence_by_nodes;
  
    weak_divergence_.resize(finite_element_space_.GetNE());
    for (int element=0; element<finite_element_space_.GetNE(); element++)
    {
      int num_dof = finite_element_space_.GetFE(element)->GetDof();
      weak_divergence_by_nodes.SetSize(num_dof, num_dof*dim_);
      weak_divergence_integrator.AssembleElementMatrix2(*finite_element_space_.GetFE(element), *finite_element_space_.GetFE(element),
                                                        *finite_element_space_.GetElementTransformation(element),
                                                        weak_divergence_by_nodes);
      weak_divergence_[element].SetSize(num_dof, num_dof*dim_);

      // Reorder so that trial space is ByDim.
      // This makes applying weak divergence to flux value simpler.
      for (int jdof=0; jdof<num_dof; jdof++)
      {
        for (int d=0; d<dim_; d++)
        {
          weak_divergence_[element].SetCol(jdof*dim_ + d, weak_divergence_by_nodes.GetColumn(d*num_dof + jdof));
        }
      }
    }
  }

  void DGAssembly::applyInverseMass(const mfem::Vector &values, mfem::Vector &dofs) const {

    for (int element=0; element<finite_element_space_.GetNE(); element++) {
      const int num_dof = finite_element_space_.GetFE(element)->GetDof();

      mfem::Array<int> vector_dofs;
      finite_element_space_.GetElementVDofs(element, vector_dofs);

      mfem::Vector element_values;
      values.GetSubVector(vector_dofs, element_values);
      mfem::DenseMatrix element_values_as_matrix(element_values.GetData(), num_dof, num_equations_);
  
      mfem::DenseMatrix result;
      result.SetSize(num_dof, num_equations_);
      mfem::Mult(inverse_mass_[element], element_values_as_matrix, result);
      dofs.SetSubVector(vector_dofs, result.GetData());
    }
  }

  void DGAssembly::computeFluxVectorInElement(const int element, const mfem::DenseMatrix &state_in_current_element, mfem::DenseMatrix &flux_in_current_element) const
  {
    mfem::Vector state_at_current_dof;
    mfem::DenseMatrix flux_at_current_dof;
    mfem::ElementTransformation* current_element_transformation = finite_element_space_.GetElementTransformation(element);
    int num_dof = finite_element_space_.GetFE(element)->GetDof();
    for (int jdof=0; jdof<num_dof; jdof++)
    {
      state_in_current_element.GetRow(jdof, state_at_current_dof);
      flux_at_current_dof.UseExternalData(flux_in_current_element.GetData() + num_equations_*dim_*jdof,
                                          num_equations_, dim_);
      flux_function_->ComputeFlux(state_at_current_dof, *current_element_transformation, flux_at_current_dof);
    }
  }

  void DGAssembly::applyWeakDivergenceInElement(
    const int element,
    const mfem::DenseMatrix &flux_in_current_element,
    mfem::DenseMatrix &rhs_in_current_element) const
  {
    mfem::AddMult_a_ABt(1.0, weak_divergence_[element], flux_in_current_element, rhs_in_current_element);
  }

  void DGAssembly::computeFluxOnElementBoundaries(const mfem::Vector &dofs, mfem::Vector &rhs) const
  {
    nonlinear_form_->AddMult(dofs, rhs);
  }

  void DGAssembly::computeHyperbolicFluxes(const mfem::Vector &dofs, mfem::Vector &rhs) const
  {
    hyperbolic_form_integrator_->ResetMaxCharSpeed();

    this->computeFluxOnElementBoundaries(dofs, rhs);

    mfem::Array<int> vector_dofs;
    mfem::DenseMatrix flux_in_current_element;
    mfem::DenseMatrix state_in_current_element;
    mfem::Vector vector_state_in_current_element;

    mfem::DenseMatrix rhs_in_current_element;
    mfem::Vector vector_rhs_in_current_element;

    for (int element=0; element<finite_element_space_.GetNE(); element++) {
      int num_dof = finite_element_space_.GetFE(element)->GetDof();
      finite_element_space_.GetElementVDofs(element, vector_dofs);

      dofs.GetSubVector(vector_dofs, vector_state_in_current_element);
      state_in_current_element.UseExternalData(vector_state_in_current_element.GetData(), num_dof, num_equations_);
      flux_in_current_element.SetSize(num_equations_, dim_*num_dof);
      this->computeFluxVectorInElement(element,state_in_current_element,flux_in_current_element);

      rhs.GetSubVector(vector_dofs, vector_rhs_in_current_element);
      rhs_in_current_element.UseExternalData(vector_rhs_in_current_element.GetData(), num_dof, num_equations_);
      flux_in_current_element.SetSize(num_equations_, dim_*num_dof);
      this->applyWeakDivergenceInElement(element, flux_in_current_element, rhs_in_current_element);

      rhs.SetSubVector(vector_dofs, vector_rhs_in_current_element);
    }
  }

  DGAssembly::~DGAssembly() = default;

} // namespace mfpic
