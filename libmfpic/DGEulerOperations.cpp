#include <libmfpic/Constants.hpp>
#include <libmfpic/DGAssembly.hpp>
#include <libmfpic/DGEulerAssembly.hpp>
#include <libmfpic/DGEulerOperations.hpp>
#include <libmfpic/ElectromagneticFieldsEvaluator.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/LowFidelityState.hpp>
#include <libmfpic/MFEMHelpers.hpp>
#include <libmfpic/Species.hpp>

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
  const LowFidelityState& state,
  const ElectromagneticFieldsEvaluator& field_evaluator) const
{
  LowFidelityState updated_state(state);

  for (int ispecies = 0; ispecies < state.numSpecies(); ++ispecies) {
    const LowFidelitySpeciesState& species_state = state.getSpeciesState(ispecies);
    LowFidelitySpeciesState& updated_species_state = updated_state.getSpeciesState(ispecies);

    temp_vector_ = 0.;
    constexpr bool include_energy_source = false;
    dg_assemblers_[ispecies]->computeElectromagneticSources(species_state, field_evaluator, temp_vector_, include_energy_source);
    dg_assemblers_[ispecies]->applyInverseMass(temp_vector_, rhs_);

    temp_vector_ = 0.;
    dg_assemblers_[ispecies]->computeIntegratedKineticEnergy(species_state, temp_vector_);
    temp_vector_ *= -1.;

    // TODO BWR at some point we may not want to deep copy updated_state so this is more safe

    mfem::GridFunction& updated_species_grid_function = updated_species_state.getGridFunction();
    updated_species_grid_function.Add(dt, rhs_);

    dg_assemblers_[ispecies]->computeIntegratedKineticEnergy(updated_species_state, temp_vector_);
    dg_assemblers_[ispecies]->applyInverseMass(temp_vector_, rhs_);

    updated_species_grid_function.Add(1., rhs_);
  }

  return updated_state;
}

LowFidelityState DGEulerOperations::move(double dt, const LowFidelityState& state) const {
  LowFidelityState updated_state(state);

  for (int ispecies = 0; ispecies < state.numSpecies(); ++ispecies) {
    const LowFidelitySpeciesState& species_state = state.getSpeciesState(ispecies);
    LowFidelitySpeciesState& updated_species_state = updated_state.getSpeciesState(ispecies);

    const mfem::GridFunction& species_grid_function = species_state.getGridFunction();

    temp_vector_ = 0.;
    dg_assemblers_[ispecies]->computeHyperbolicFluxes(species_grid_function, temp_vector_);
    dg_assemblers_[ispecies]->applyInverseMass(temp_vector_, rhs_);

    mfem::GridFunction& updated_species_grid_function = updated_species_state.getGridFunction();
    updated_species_grid_function.Add(dt, rhs_);
  }

  return updated_state;
}

LowFidelityState DGEulerOperations::moveAccelerate(
  const double dt,
  const LowFidelityState& state,
  const ElectromagneticFieldsEvaluator& field_evaluator) const
{
  LowFidelityState updated_state(state);
  for (int i_species = 0; i_species < state.numSpecies(); ++i_species) {
    const LowFidelitySpeciesState& species_state = state.getSpeciesState(i_species);
    const mfem::GridFunction& species_grid_function = species_state.getGridFunction();

    temp_vector_ = 0.;
    dg_assemblers_[i_species]->computeHyperbolicFluxes(species_grid_function, temp_vector_);
    constexpr bool include_energy_source = true;
    dg_assemblers_[i_species]->computeElectromagneticSources(species_state, field_evaluator, temp_vector_, include_energy_source);
    dg_assemblers_[i_species]->applyInverseMass(temp_vector_, rhs_);

    LowFidelitySpeciesState& updated_species_state = updated_state.getSpeciesState(i_species);
    mfem::GridFunction& updated_species_grid_function = updated_species_state.getGridFunction();
    updated_species_grid_function.Add(dt, rhs_);
  }
  return updated_state;
}

LowFidelityState DGEulerOperations::computeRHS(
  const LowFidelityState& state,
  const ElectromagneticFieldsEvaluator& field_evaluator) const
{
  LowFidelityState rhs(state);
  for (int i_species = 0; i_species < state.numSpecies(); ++i_species) {
    const LowFidelitySpeciesState& species_state = state.getSpeciesState(i_species);
    const mfem::GridFunction& species_grid_function = species_state.getGridFunction();

    temp_vector_ = 0.;
    dg_assemblers_[i_species]->computeHyperbolicFluxes(species_grid_function, temp_vector_);
    constexpr bool include_energy_source = true;
    dg_assemblers_[i_species]->computeElectromagneticSources(species_state, field_evaluator, temp_vector_, include_energy_source);
    dg_assemblers_[i_species]->applyInverseMass(temp_vector_, rhs_);

    LowFidelitySpeciesState& rhs_species_state = rhs.getSpeciesState(i_species);
    mfem::GridFunction& rhs_species_grid_function = rhs_species_state.getGridFunction();
    rhs_species_grid_function.Set(1., rhs_);
  }
  return rhs;
}

  LowFidelityState DGEulerOperations::addVolumetricSource(
    double dt,
    const LowFidelityState& state) const
  {
    LowFidelityState updated_state(state);

    for (int ispecies = 0; ispecies < state.numSpecies(); ++ispecies) {
      const LowFidelitySpeciesState& species_state = state.getSpeciesState(ispecies);
      LowFidelitySpeciesState& updated_species_state = updated_state.getSpeciesState(ispecies);

      temp_vector_ = 0.;
      dg_assemblers_[ispecies]->computeVolumetricSources(species_state, temp_vector_);
      dg_assemblers_[ispecies]->applyInverseMass(temp_vector_, rhs_);

      mfem::GridFunction& updated_species_grid_function = updated_species_state.getGridFunction();
      updated_species_grid_function.Add(dt, rhs_);
    }

    return updated_state;
  }

IntegratedCharge DGEulerOperations::assembleCharge(const LowFidelityState& state) const {
  IntegratedCharge charge_state(charge_discretization_);
  charge_state.setIntegratedChargeValue(0.0);

  mfem::FiniteElementSpace & finite_element_space = charge_discretization_.getFeSpace();

  for (int ispecies = 0; ispecies < state.numSpecies(); ++ispecies) {
    const LowFidelitySpeciesState& species_state = state.getSpeciesState(ispecies);

    Species species = species_state.getSpecies();
    // don't use species.charge_over_mass in order to support ghost species
    double charge_over_mass = species.charge / species.mass;

    const mfem::GridFunction& species_grid_function = species_state.getGridFunction();

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
      species_grid_function.GetVectorValues(
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
            weight * basis_values(jdof) * charge_over_mass * fluid_state(euler::ConservativeVariables::MASS_DENSITY));
        }
      }

    }
  }
  return charge_state;
}

double DGEulerOperations::estimateCFL(const double & dt, const double & smallest_cell_lengthscale) const {
  double max_speed = 0.;
  for (size_t ispecies = 0; ispecies < dg_assemblers_.size(); ++ispecies) {
    max_speed = fmax(max_speed, dg_assemblers_[ispecies]->getMaxCharSpeed());
  }

  return max_speed * dt / smallest_cell_lengthscale;
}

double DGEulerOperations::computeTotalEnergy(const LowFidelityState& state) const {
  Discretization identity_test_function_discretization = getIdentityTestFunctionDiscretization(charge_discretization_);

  double total_energy = 0;
  for (int i_species = 0; i_species < state.numSpecies(); ++i_species) {
    const LowFidelitySpeciesState& species_state = state.getSpeciesState(i_species);
    const mfem::GridFunction& grid_function = species_state.getGridFunction();

    // GridFunctionCoefficient component is 1 based not 0 based indexing
    constexpr int component = euler::ConservativeVariables::TOTAL_ENERGY_DENSITY + 1;
    mfem::GridFunctionCoefficient species_total_energy_density_coefficient(&grid_function, component);

    mfem::LinearForm species_total_energy_by_cell(&identity_test_function_discretization.getFeSpace());
    species_total_energy_by_cell.AddDomainIntegrator(new mfem::DomainLFIntegrator(species_total_energy_density_coefficient));
    species_total_energy_by_cell.Assemble();
    total_energy += species_total_energy_by_cell.Sum();
  }

  return total_energy;
}

double DGEulerOperations::computeTotalKineticEnergy(const LowFidelityState& state) const {
  Discretization identity_test_function_discretization = getIdentityTestFunctionDiscretization(charge_discretization_);

  double total_kinetic_energy = 0.;
  for (int i_species = 0; i_species < state.numSpecies(); ++i_species) {
    const LowFidelitySpeciesState& species_state = state.getSpeciesState(i_species);
    const mfem::GridFunction& grid_function = species_state.getGridFunction();

    auto conservative_state_coefficient = std::make_unique<mfem::VectorGridFunctionCoefficient>(&grid_function);
    TransformedVectorCoefficient kinetic_energy_density_coefficient(
      std::move(conservative_state_coefficient),
      euler::getKineticEnergyDensityFromConservativeState);

    mfem::LinearForm species_kinetic_energy_by_cell(&identity_test_function_discretization.getFeSpace());
    species_kinetic_energy_by_cell.AddDomainIntegrator(new mfem::DomainLFIntegrator(kinetic_energy_density_coefficient));
    species_kinetic_energy_by_cell.Assemble();

    const double species_total_kinetic_energy = species_kinetic_energy_by_cell.Sum();
    total_kinetic_energy += species_total_kinetic_energy;
  }

  return total_kinetic_energy;
}

double DGEulerOperations::computeTotalCharge(const LowFidelityState& state) const {
  IntegratedCharge integrated_charge = this->assembleCharge(state);
  mfem::Vector integrated_charge_vector = integrated_charge.getIntegratedCharge();
  double sum = 0.;
  for (int i = 0; i < integrated_charge_vector.Size(); ++i) {
    sum += integrated_charge_vector(i);
  }
  return sum;
}

double DGEulerOperations::evaluateParticleDistributionFunction(
  const LowFidelityState& state,
  const mfem::Vector position,
  const mfem::Vector velocity,
  const int element,
  const Species& species_to_evaluate) const
{
  mfem::FiniteElementSpace & finite_element_space = charge_discretization_.getFeSpace();
  mfem::Mesh * mesh = finite_element_space.GetMesh();
  for (int ispecies = 0; ispecies < state.numSpecies(); ++ispecies) {
    const LowFidelitySpeciesState& species_state = state.getSpeciesState(ispecies);
    Species species = species_state.getSpecies();
    if (species == species_to_evaluate) {
      const mfem::GridFunction& species_grid_function = species_state.getGridFunction();
      mfem::ElementTransformation *element_transformation = mesh->GetElementTransformation(element);

      mfem::InverseElementTransformation inverse_element_transformation(element_transformation);
      mfem::IntegrationPoint ip_ref;
      int info = inverse_element_transformation.Transform(position, ip_ref);
      MFEM_VERIFY(info == mfem::InverseElementTransformation::Inside,
          "Point is not inside the element.");
      mfem::Vector fluid_state_at_position;
      species_grid_function.GetVectorValue(element, ip_ref, fluid_state_at_position);

      mfem::Vector primitive_state = euler::convertFromConservativeToPrimitive(fluid_state_at_position, species);
      double particle_distribution_function_value = euler::evaluateMaxwellian(primitive_state,velocity, species);
      return particle_distribution_function_value;
    }
  }
  std::ostringstream error_message;
  error_message << "Species not found in low fidelity state.\n";
  errorWithUserMessage(error_message.str()); 
}

std::unordered_map<Species,mfem::Vector> DGEulerOperations::integralForVarianceReducedNumberDensity(mfem::FiniteElementSpace finite_element_space, const LowFidelityState& current_state) const
  {
    mfem::Mesh& mesh = *finite_element_space.GetMesh();
    std::unordered_map<Species, mfem::Vector> number_density_integral;

    for (int ispecies = 0; ispecies < current_state.numSpecies(); ++ispecies) {
      const LowFidelitySpeciesState& current_species_state = current_state.getSpeciesState(ispecies);
      Species current_species = current_species_state.getSpecies();
      number_density_integral.insert({current_species, mfem::Vector(mesh.GetNE())});
      number_density_integral.at(current_species) = 0.0;
      const mfem::GridFunction& current_species_grid_function = current_species_state.getGridFunction();
      mfem::DenseMatrix fluid_state_at_integration_point_locations, integration_point_locations_in_physical_frame;

      for (int element=0; element<finite_element_space.GetNE(); element++)
      {
        const mfem::IntegrationRule &integration_rule = mfem::IntRules.Get(
          finite_element_space.GetFE(element)->GetGeomType(),
          2*finite_element_space.GetFE(element)->GetOrder());

        mfem::ElementTransformation* element_transformation = finite_element_space.GetElementTransformation(element);
        element_transformation->Transform(integration_rule, integration_point_locations_in_physical_frame);
        current_species_grid_function.GetVectorValues(
          *element_transformation, integration_rule, fluid_state_at_integration_point_locations);

        mfem::Vector position(integration_point_locations_in_physical_frame.NumRows());
        mfem::Vector fluid_state(dg_assemblers_[ispecies]->getNumberOfEquations());

        for (int ipoint = 0; ipoint < integration_rule.GetNPoints(); ++ipoint) 
        {
          const mfem::IntegrationPoint &integration_point = integration_rule.IntPoint(ipoint); 
          element_transformation->SetIntPoint(&integration_point);
          integration_point_locations_in_physical_frame.GetColumn(ipoint, position);
          fluid_state_at_integration_point_locations.GetColumn(ipoint, fluid_state);
          mfem::Vector primitive_state = euler::convertFromConservativeToPrimitive(fluid_state, current_species);
          const double weight = integration_point.weight * element_transformation->Weight();
          number_density_integral.at(current_species)(element) += weight * primitive_state(euler::PrimitiveVariables::NUMBER_DENSITY);
        }
      }
    }
    return number_density_integral;
  }

  std::unordered_map<Species,mfem::DenseMatrix> DGEulerOperations::integralForVarianceReducedBulkVelocity(mfem::FiniteElementSpace finite_element_space, const LowFidelityState& current_state) const
  {
    mfem::Mesh& mesh = *finite_element_space.GetMesh();
    std::unordered_map<Species, mfem::DenseMatrix> bulk_velocity_integral;

    for (int ispecies = 0; ispecies < current_state.numSpecies(); ++ispecies) {
      const LowFidelitySpeciesState& current_species_state = current_state.getSpeciesState(ispecies);
      Species current_species = current_species_state.getSpecies();
      bulk_velocity_integral.insert({current_species, mfem::DenseMatrix(3,mesh.GetNE())});
      bulk_velocity_integral.at(current_species) = 0.0;
      const mfem::GridFunction& current_species_grid_function = current_species_state.getGridFunction();
      mfem::DenseMatrix fluid_state_at_integration_point_locations, integration_point_locations_in_physical_frame;

      for (int element=0; element<finite_element_space.GetNE(); element++)
      {
        const mfem::IntegrationRule &integration_rule = mfem::IntRules.Get(
          finite_element_space.GetFE(element)->GetGeomType(),
          2*finite_element_space.GetFE(element)->GetOrder());

        mfem::ElementTransformation* element_transformation = finite_element_space.GetElementTransformation(element);
        element_transformation->Transform(integration_rule, integration_point_locations_in_physical_frame);
        current_species_grid_function.GetVectorValues(
          *element_transformation, integration_rule, fluid_state_at_integration_point_locations);

        mfem::Vector position(integration_point_locations_in_physical_frame.NumRows());
        mfem::Vector fluid_state(dg_assemblers_[ispecies]->getNumberOfEquations());

        for (int ipoint = 0; ipoint < integration_rule.GetNPoints(); ++ipoint) 
        {
          const mfem::IntegrationPoint &integration_point = integration_rule.IntPoint(ipoint); 
          element_transformation->SetIntPoint(&integration_point);
          integration_point_locations_in_physical_frame.GetColumn(ipoint, position);
          fluid_state_at_integration_point_locations.GetColumn(ipoint, fluid_state);
          mfem::Vector primitive_state = euler::convertFromConservativeToPrimitive(fluid_state, current_species);
          const double weight = integration_point.weight * element_transformation->Weight();
          bulk_velocity_integral.at(current_species)(0,element) += weight * primitive_state(euler::PrimitiveVariables::NUMBER_DENSITY) * primitive_state(euler::PrimitiveVariables::X_BULK_VELOCITY);
          bulk_velocity_integral.at(current_species)(1,element) += weight * primitive_state(euler::PrimitiveVariables::NUMBER_DENSITY) * primitive_state(euler::PrimitiveVariables::Y_BULK_VELOCITY);
          bulk_velocity_integral.at(current_species)(2,element) += weight * primitive_state(euler::PrimitiveVariables::NUMBER_DENSITY) * primitive_state(euler::PrimitiveVariables::Z_BULK_VELOCITY);
        }
      }
    }
    return bulk_velocity_integral;
  }

  std::unordered_map<Species,mfem::Vector> DGEulerOperations::integralForVarianceReducedTemperature(mfem::FiniteElementSpace finite_element_space, const LowFidelityState& current_state) const
  {
    mfem::Mesh& mesh = *finite_element_space.GetMesh();
    std::unordered_map<Species, mfem::Vector> temperature_integral;

    for (int ispecies = 0; ispecies < current_state.numSpecies(); ++ispecies) {
      const LowFidelitySpeciesState& current_species_state = current_state.getSpeciesState(ispecies);
      Species current_species = current_species_state.getSpecies();
      temperature_integral.insert({current_species, mfem::Vector(mesh.GetNE())});
      temperature_integral.at(current_species) = 0.0;
      const mfem::GridFunction& current_species_grid_function = current_species_state.getGridFunction();
      mfem::DenseMatrix fluid_state_at_integration_point_locations, integration_point_locations_in_physical_frame;

      for (int element=0; element<finite_element_space.GetNE(); element++)
      {
        const mfem::IntegrationRule &integration_rule = mfem::IntRules.Get(
          finite_element_space.GetFE(element)->GetGeomType(),
          2*finite_element_space.GetFE(element)->GetOrder());

        mfem::ElementTransformation* element_transformation = finite_element_space.GetElementTransformation(element);
        element_transformation->Transform(integration_rule, integration_point_locations_in_physical_frame);
        current_species_grid_function.GetVectorValues(
          *element_transformation, integration_rule, fluid_state_at_integration_point_locations);

        mfem::Vector position(integration_point_locations_in_physical_frame.NumRows());
        mfem::Vector fluid_state(dg_assemblers_[ispecies]->getNumberOfEquations());

        for (int ipoint = 0; ipoint < integration_rule.GetNPoints(); ++ipoint) 
        {
          const mfem::IntegrationPoint &integration_point = integration_rule.IntPoint(ipoint); 
          element_transformation->SetIntPoint(&integration_point);
          integration_point_locations_in_physical_frame.GetColumn(ipoint, position);
          fluid_state_at_integration_point_locations.GetColumn(ipoint, fluid_state);
          mfem::Vector primitive_state = euler::convertFromConservativeToPrimitive(fluid_state, current_species);
          const double weight = integration_point.weight * element_transformation->Weight();
          const double bulk_velocity_mag_squared
            = primitive_state(euler::PrimitiveVariables::X_BULK_VELOCITY) *primitive_state(euler::PrimitiveVariables::X_BULK_VELOCITY)
            + primitive_state(euler::PrimitiveVariables::Y_BULK_VELOCITY) *primitive_state(euler::PrimitiveVariables::Y_BULK_VELOCITY)
            + primitive_state(euler::PrimitiveVariables::Z_BULK_VELOCITY) *primitive_state(euler::PrimitiveVariables::Z_BULK_VELOCITY);
          temperature_integral.at(current_species)(element) += weight * primitive_state(euler::PrimitiveVariables::NUMBER_DENSITY) * (primitive_state(euler::PrimitiveVariables::TEMPERATURE) + current_species.mass/(3 * constants::boltzmann_constant) * bulk_velocity_mag_squared);
        }
      }
    }
    return temperature_integral;
  }


std::unordered_map<Species,mfem::Vector> DGEulerOperations::getCellAveragedNumberDensity(const LowFidelityState& current_state) const
{
  mfem::FiniteElementSpace finite_element_space = charge_discretization_.getFeSpace();
  mfem::Mesh& mesh = *finite_element_space.GetMesh();
  std::unordered_map<Species,mfem::Vector> number_density_integral;

  for (int ispecies = 0; ispecies < current_state.numSpecies(); ++ispecies) {
    const LowFidelitySpeciesState& current_species_state = current_state.getSpeciesState(ispecies);
    Species current_species = current_species_state.getSpecies();
    number_density_integral.insert({current_species, mfem::Vector(mesh.GetNE())});
    number_density_integral.at(current_species) = 0.0;
    const mfem::GridFunction& current_species_grid_function = current_species_state.getGridFunction();
    mfem::DenseMatrix fluid_state_at_integration_point_locations, integration_point_locations_in_physical_frame;

    for (int element=0; element<finite_element_space.GetNE(); element++)
    {
      const mfem::IntegrationRule &integration_rule = mfem::IntRules.Get(
        finite_element_space.GetFE(element)->GetGeomType(),
        2*finite_element_space.GetFE(element)->GetOrder());

      mfem::ElementTransformation* element_transformation = finite_element_space.GetElementTransformation(element);
      element_transformation->Transform(integration_rule, integration_point_locations_in_physical_frame);
      current_species_grid_function.GetVectorValues(
        *element_transformation, integration_rule, fluid_state_at_integration_point_locations);

      mfem::Vector position(integration_point_locations_in_physical_frame.NumRows());
      mfem::Vector fluid_state(dg_assemblers_[ispecies]->getNumberOfEquations());
      const double element_volume = mesh.GetElementVolume(element);

      for (int ipoint = 0; ipoint < integration_rule.GetNPoints(); ++ipoint) 
      {
        const mfem::IntegrationPoint &integration_point = integration_rule.IntPoint(ipoint); 
        element_transformation->SetIntPoint(&integration_point);
        integration_point_locations_in_physical_frame.GetColumn(ipoint, position);
        fluid_state_at_integration_point_locations.GetColumn(ipoint, fluid_state);
        mfem::Vector primitive_state = euler::convertFromConservativeToPrimitive(fluid_state, current_species);
        const double weight = integration_point.weight * element_transformation->Weight();
        number_density_integral.at(current_species)(element) += weight * primitive_state(euler::PrimitiveVariables::NUMBER_DENSITY) / element_volume;
      }
    }
  }
  return number_density_integral;
}

std::unordered_map<Species,mfem::DenseMatrix> DGEulerOperations::getCellAveragedBulkVelocity(const LowFidelityState& current_state) const
{
  mfem::FiniteElementSpace finite_element_space = charge_discretization_.getFeSpace();
  mfem::Mesh& mesh = *finite_element_space.GetMesh();
  std::unordered_map<Species, mfem::DenseMatrix> bulk_velocity_integral;

  for (int ispecies = 0; ispecies < current_state.numSpecies(); ++ispecies) {
    const LowFidelitySpeciesState& current_species_state = current_state.getSpeciesState(ispecies);
    Species current_species = current_species_state.getSpecies();
    bulk_velocity_integral.insert({current_species, mfem::DenseMatrix(3,mesh.GetNE())});
    bulk_velocity_integral.at(current_species) = 0.0;
    const mfem::GridFunction& current_species_grid_function = current_species_state.getGridFunction();
    mfem::DenseMatrix fluid_state_at_integration_point_locations, integration_point_locations_in_physical_frame;

    for (int element=0; element<finite_element_space.GetNE(); element++)
    {
      const mfem::IntegrationRule &integration_rule = mfem::IntRules.Get(
        finite_element_space.GetFE(element)->GetGeomType(),
        2*finite_element_space.GetFE(element)->GetOrder());

      mfem::ElementTransformation* element_transformation = finite_element_space.GetElementTransformation(element);
      element_transformation->Transform(integration_rule, integration_point_locations_in_physical_frame);
      current_species_grid_function.GetVectorValues(
        *element_transformation, integration_rule, fluid_state_at_integration_point_locations);

      mfem::Vector position(integration_point_locations_in_physical_frame.NumRows());
      mfem::Vector fluid_state(dg_assemblers_[ispecies]->getNumberOfEquations());
      const double element_volume = mesh.GetElementVolume(element);

      for (int ipoint = 0; ipoint < integration_rule.GetNPoints(); ++ipoint) 
      {
        const mfem::IntegrationPoint &integration_point = integration_rule.IntPoint(ipoint); 
        element_transformation->SetIntPoint(&integration_point);
        integration_point_locations_in_physical_frame.GetColumn(ipoint, position);
        fluid_state_at_integration_point_locations.GetColumn(ipoint, fluid_state);
        mfem::Vector primitive_state = euler::convertFromConservativeToPrimitive(fluid_state, current_species);
        const double weight = integration_point.weight * element_transformation->Weight();
        mfem::Vector velocity_in_element(bulk_velocity_integral.at(current_species).GetColumn(element), 3);
        velocity_in_element(0) += weight * primitive_state(euler::PrimitiveVariables::X_BULK_VELOCITY) / element_volume;
        velocity_in_element(1) += weight * primitive_state(euler::PrimitiveVariables::Y_BULK_VELOCITY) / element_volume;
        velocity_in_element(2) += weight * primitive_state(euler::PrimitiveVariables::Z_BULK_VELOCITY) / element_volume;
      }
    }
  }
  return bulk_velocity_integral;
}

std::unordered_map<Species,mfem::Vector> DGEulerOperations::getCellAveragedTemperature(const LowFidelityState& current_state) const
{
  mfem::FiniteElementSpace finite_element_space = charge_discretization_.getFeSpace();
  mfem::Mesh& mesh = *finite_element_space.GetMesh();
  std::unordered_map<Species, mfem::Vector> temperature_integral;

  for (int ispecies = 0; ispecies < current_state.numSpecies(); ++ispecies) {
    const LowFidelitySpeciesState& current_species_state = current_state.getSpeciesState(ispecies);
    Species current_species = current_species_state.getSpecies();
    temperature_integral.insert({current_species, mfem::Vector(mesh.GetNE())});
    temperature_integral.at(current_species) = 0.0;
    const mfem::GridFunction& current_species_grid_function = current_species_state.getGridFunction();
    mfem::DenseMatrix fluid_state_at_integration_point_locations, integration_point_locations_in_physical_frame;

    for (int element=0; element<finite_element_space.GetNE(); element++)
    {
      const mfem::IntegrationRule &integration_rule = mfem::IntRules.Get(
        finite_element_space.GetFE(element)->GetGeomType(),
        2*finite_element_space.GetFE(element)->GetOrder());

      mfem::ElementTransformation* element_transformation = finite_element_space.GetElementTransformation(element);
      element_transformation->Transform(integration_rule, integration_point_locations_in_physical_frame);
      current_species_grid_function.GetVectorValues(
        *element_transformation, integration_rule, fluid_state_at_integration_point_locations);

      mfem::Vector position(integration_point_locations_in_physical_frame.NumRows());
      mfem::Vector fluid_state(dg_assemblers_[ispecies]->getNumberOfEquations());
      const double element_volume = mesh.GetElementVolume(element);

      for (int ipoint = 0; ipoint < integration_rule.GetNPoints(); ++ipoint) 
      {
        const mfem::IntegrationPoint &integration_point = integration_rule.IntPoint(ipoint); 
        element_transformation->SetIntPoint(&integration_point);
        integration_point_locations_in_physical_frame.GetColumn(ipoint, position);
        fluid_state_at_integration_point_locations.GetColumn(ipoint, fluid_state);
        mfem::Vector primitive_state = euler::convertFromConservativeToPrimitive(fluid_state, current_species);
        const double weight = integration_point.weight * element_transformation->Weight();
        temperature_integral.at(current_species)(element) += weight * primitive_state(euler::PrimitiveVariables::TEMPERATURE) / element_volume;
      }
    }
  }
  return temperature_integral;
}

} // namespace mfpic
