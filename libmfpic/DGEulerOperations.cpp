#include <libmfpic/DGAssembly.hpp>
#include <libmfpic/DGEulerAssembly.hpp>
#include <libmfpic/DGEulerOperations.hpp>
#include <libmfpic/ElectromagneticFieldsEvaluator.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/LowFidelityState.hpp>

namespace mfpic {

  DGEulerOperations::DGEulerOperations(
    Discretization &charge_discretization,
    std::vector<std::shared_ptr<DGEulerAssembly>> &dg_assemblers
    ) :
    charge_discretization_(charge_discretization),
    dg_assemblers_(dg_assemblers)
  { 
    rhs_ = dg_assemblers[0]->getFEVector();
    temp_vector_ = dg_assemblers[0]->getFEVector(); 
  }

  LowFidelityState DGEulerOperations::accelerate(
    double dt,
    const LowFidelityState& current_state,
    const ElectromagneticFieldsEvaluator& field_evaluator) const
  {
    auto& message = std::cout;
    message << "DGEulerOperations::accelerate" << std::endl;
    message << "current_state.numSpecies() = " << current_state.numSpecies() << std::endl;

    LowFidelityState updated_state(current_state);

    for (int ispecies = 0; ispecies < current_state.numSpecies(); ++ispecies) {
      message << "ispecies = " << ispecies << std::endl;
      const LowFidelitySpeciesState& current_species_state = current_state.getSpeciesState(ispecies);
      LowFidelitySpeciesState& updated_species_state = updated_state.getSpeciesState(ispecies);

      message << "compute sources" << std::endl;
      temp_vector_ = 0.;
      dg_assemblers_[ispecies]->computeSources(current_species_state, field_evaluator, temp_vector_);
      message << "applyInverseMass" << std::endl;
      dg_assemblers_[ispecies]->applyInverseMass(temp_vector_, rhs_);

      temp_vector_ = 0.;
      message << "computeIntegratedKineticEnergy" << std::endl;
      dg_assemblers_[ispecies]->computeIntegratedKineticEnergy(current_species_state, temp_vector_);
      temp_vector_ *= -1.;

      // TODO BWR at some point we may not want to deep copy updated_state so this is more safe

      message << "add" << std::endl;
      mfem::GridFunction& updated_species_grid_function = updated_species_state.getGridFunction();
      updated_species_grid_function.Add(dt, rhs_);

      message << "computeIntegratedKineticEnergy" << std::endl;
      dg_assemblers_[ispecies]->computeIntegratedKineticEnergy(updated_species_state, temp_vector_);
      message << "applyInverseMass" << std::endl;
      dg_assemblers_[ispecies]->applyInverseMass(temp_vector_, rhs_);

      message << "Add" << std::endl;
      updated_species_grid_function.Add(1., rhs_);
    }

    return updated_state;
  }

  LowFidelityState DGEulerOperations::move(double dt, const LowFidelityState& current_state) const
  {
    LowFidelityState updated_state(current_state);

    for (int ispecies = 0; ispecies < current_state.numSpecies(); ++ispecies) {
      const LowFidelitySpeciesState& current_species_state = current_state.getSpeciesState(ispecies);
      LowFidelitySpeciesState& updated_species_state = updated_state.getSpeciesState(ispecies);

      const mfem::GridFunction& current_species_grid_function = current_species_state.getGridFunction();

      temp_vector_ = 0.;
      dg_assemblers_[ispecies]->computeHyperbolicFluxes(current_species_grid_function, temp_vector_);
      dg_assemblers_[ispecies]->applyInverseMass(temp_vector_, rhs_);

      mfem::GridFunction& updated_species_grid_function = updated_species_state.getGridFunction();
      updated_species_grid_function.Add(dt, rhs_);
    }

    return updated_state;
  }

  IntegratedCharge DGEulerOperations::assembleCharge(const LowFidelityState& current_state) const
  {
    IntegratedCharge charge_state(charge_discretization_);
    charge_state.setIntegratedChargeValue(0.0);

    mfem::FiniteElementSpace & finite_element_space = charge_discretization_.getFeSpace();

    for (int ispecies = 0; ispecies < current_state.numSpecies(); ++ispecies) {
      const LowFidelitySpeciesState& current_species_state = current_state.getSpeciesState(ispecies);

      Species current_species = current_species_state.getSpecies();
      double species_charge_over_mass = current_species.charge / current_species.mass;

      const mfem::GridFunction& current_species_grid_function = current_species_state.getGridFunction();

      mfem::DenseMatrix fluid_state_at_integration_point_locations, integration_point_locations_in_physical_frame;
      mfem::Array<int> vector_dofs;

      for (int element=0; element<finite_element_space.GetNE(); element++)
      {
        finite_element_space.GetElementVDofs(element, vector_dofs);

        const mfem::IntegrationRule &integration_rule = mfem::IntRules.Get(
          finite_element_space.GetFE(element)->GetGeomType(),
          2*finite_element_space.GetFE(element)->GetOrder());

        mfem::ElementTransformation* element_transformation = finite_element_space.GetElementTransformation(element);
        int num_element_dof = finite_element_space.GetFE(element)->GetDof();
        element_transformation->Transform(integration_rule, integration_point_locations_in_physical_frame);
        current_species_grid_function.GetVectorValues(
          *element_transformation, integration_rule, fluid_state_at_integration_point_locations);

        mfem::Vector basis_values(num_element_dof);
        mfem::Vector position(integration_point_locations_in_physical_frame.NumRows());
        mfem::Vector fluid_state(dg_assemblers_[ispecies]->getNumberOfEquations());

        for (int ipoint = 0; ipoint < integration_rule.GetNPoints(); ++ipoint) 
        {
          const mfem::IntegrationPoint &integration_point = integration_rule.IntPoint(ipoint); 
          element_transformation->SetIntPoint(&integration_point);
          finite_element_space.GetFE(element)->CalcShape(integration_point, basis_values);
          integration_point_locations_in_physical_frame.GetColumn(ipoint, position);
          fluid_state_at_integration_point_locations.GetColumn(ipoint, fluid_state);
          const double weight = integration_point.weight * element_transformation->Weight();
          for (int jdof=0; jdof<num_element_dof; jdof++)
          {
            charge_state.addIntegratedChargeValue(
              vector_dofs[jdof],
              weight * basis_values(jdof) * species_charge_over_mass * fluid_state(euler::ConservativeVariables::MASS_DENSITY));
          }
        }

      }
    }
    return charge_state;
  }

  double DGEulerOperations::estimateCFL(const double & dt, const double & smallest_cell_lengthscale) const 
  {

    double max_speed = 0.;
    for (size_t ispecies = 0; ispecies < dg_assemblers_.size(); ++ispecies) {
      max_speed = fmax(max_speed, dg_assemblers_[ispecies]->getMaxCharSpeed());
    }

    return max_speed * dt / smallest_cell_lengthscale;
  }

} // namespace
