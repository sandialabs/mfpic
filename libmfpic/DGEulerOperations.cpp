#include <libmfpic/Constants.hpp>
#include <libmfpic/DGAssembly.hpp>
#include <libmfpic/DGEulerAssembly.hpp>
#include <libmfpic/DGEulerOperations.hpp>
#include <libmfpic/ElectromagneticFieldsEvaluator.hpp>
#include <libmfpic/ElectrostaticFieldState.hpp>
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

std::pair<LowFidelityState, ElectrostaticFieldState> DGEulerOperations::plasmaOscillate(
  const double dt,
  const LowFidelityState& state,
  const ElectrostaticFieldState& field_state) const
{
  LowFidelityState updated_state(state);
  ElectrostaticFieldState updated_field_state(field_state);

  assert(state.numSpecies() == 2);

  const LowFidelitySpeciesState& species_state_0 = state.getSpeciesState(0);
  const LowFidelitySpeciesState& species_state_1 = state.getSpeciesState(1);
  const mfem::GridFunction& euler_grid_function_0 = species_state_0.getGridFunction();
  const mfem::GridFunction& euler_grid_function_1 = species_state_1.getGridFunction();

  LowFidelitySpeciesState& updated_species_state_0 = updated_state.getSpeciesState(0);
  LowFidelitySpeciesState& updated_species_state_1 = updated_state.getSpeciesState(1);
  mfem::GridFunction& updated_euler_grid_function_0 = updated_species_state_0.getGridFunction();
  mfem::GridFunction& updated_euler_grid_function_1 = updated_species_state_1.getGridFunction();

  const mfem::FiniteElementSpace* euler_fe_space = euler_grid_function_0.FESpace();
  assert(euler_fe_space->GetMaxElementOrder() == 0);

  const mfem::GridFunction& e_field_grid_function = field_state.getEFieldGridFunction();
  mfem::GridFunction& updated_e_field_grid_function = updated_field_state.getEFieldGridFunction();

  const mfem::FiniteElementSpace* e_field_fe_space = e_field_grid_function.FESpace();
  assert(e_field_fe_space->GetMaxElementOrder() == 0);

  const Species species_0 = species_state_0.getSpecies();
  const Species species_1 = species_state_1.getSpecies();

  const double mass_0 = species_0.mass;
  const double mass_1 = species_1.mass;

  const double charge_0 = species_0.charge;
  const double charge_1 = species_1.charge;

  for (int i_element = 0; i_element < euler_fe_space->GetNE(); ++i_element) {
    mfem::Array<int> index_of_euler_dofs_on_cell;
    euler_fe_space->GetElementVDofs(i_element, index_of_euler_dofs_on_cell);

    // using the fact dg_order = 0 which implies only one dof per equation on each cell
    mfem::Vector euler_dofs_on_cell_0;
    mfem::Vector euler_dofs_on_cell_1;
    euler_grid_function_0.GetSubVector(index_of_euler_dofs_on_cell, euler_dofs_on_cell_0);
    euler_grid_function_1.GetSubVector(index_of_euler_dofs_on_cell, euler_dofs_on_cell_1);

    mfem::Array<int> index_of_e_field_dofs_on_cell;
    e_field_fe_space->GetElementVDofs(i_element, index_of_e_field_dofs_on_cell);

    mfem::Vector En;
    e_field_grid_function.GetSubVector(index_of_e_field_dofs_on_cell, En);

    // using the fact that the value of dof is the value of the equation over the cell
    const double rho_0 = euler_dofs_on_cell_0[euler::ConservativeVariables::MASS_DENSITY];
    const double rho_1 = euler_dofs_on_cell_1[euler::ConservativeVariables::MASS_DENSITY];

    const double internal_energy_density_0 = euler::getInternalEnergyDensityFromConservativeState(euler_dofs_on_cell_0);
    const double internal_energy_density_1 = euler::getInternalEnergyDensityFromConservativeState(euler_dofs_on_cell_1);

    const mfem::Vector p0n = euler::getMomentumDensityFromConservativeState(euler_dofs_on_cell_0);
    const mfem::Vector p1n = euler::getMomentumDensityFromConservativeState(euler_dofs_on_cell_1);

    const double denominator = 4. * mass_0 * mass_0 * mass_1 * mass_1 * constants::permittivity
      + dt * dt * (mass_1 * mass_1 * charge_0 * charge_0 * rho_0 + mass_0 * mass_0 * charge_1 * charge_1 * rho_1);

    const double p0np1_En_numerator = 4. * mass_0 * mass_1 * mass_1 * charge_0 * dt * constants::permittivity * rho_0;
    const double p0np1_p0n_numerator = 4. * mass_0 * mass_0 * mass_1 * mass_1 * constants::permittivity
      - dt * dt * (mass_1 * mass_1 * charge_0 * charge_0 * rho_0 - mass_0 * mass_0 * charge_1 * charge_1 * rho_1);
    const double p0np1_p1n_numerator = -2. * mass_0 * mass_1 * charge_0 * charge_1 * dt * dt * rho_0;

    const double p1np1_En_numerator = 4. * mass_0 * mass_0 * mass_1 * charge_1 * dt * constants::permittivity * rho_1;
    const double p1np1_p0n_numerator = -2. * mass_0 * mass_1 * charge_0 * charge_1 * dt * dt * rho_1;
    const double p1np1_p1n_numerator = 4. * mass_0 * mass_0 * mass_1 * mass_1 * constants::permittivity
      + dt * dt * (mass_1 * mass_1 * charge_0 * charge_0 * rho_0 - mass_0 * mass_0 * charge_1 * charge_1 * rho_1);

    const double Enp1_En_numerator = 4. * mass_0 * mass_0 * mass_1 * mass_1 * constants::permittivity
      - dt * dt * (mass_1 * mass_1 * charge_0 * charge_0 * rho_0 + mass_0 * mass_0 * charge_1 * charge_1 * rho_1);
    const double Enp1_p0n_numerator = -4. * mass_0 * mass_1 * mass_1 * charge_0 * dt;
    const double Enp1_p1n_numerator = -4. * mass_0 * mass_0 * mass_1 * charge_1 * dt;

    mfem::Vector p0np1(En.Size());
    p0np1.Set(p0np1_En_numerator / denominator, En);
    p0np1.Add(p0np1_p0n_numerator / denominator, p0n);
    p0np1.Add(p0np1_p1n_numerator / denominator, p1n);

    mfem::Vector p1np1(En.Size());
    p1np1.Set(p1np1_En_numerator / denominator, En);
    p1np1.Add(p1np1_p0n_numerator / denominator, p0n);
    p1np1.Add(p1np1_p1n_numerator / denominator, p1n);

    mfem::Vector Enp1(En.Size());
    Enp1.Set(Enp1_En_numerator / denominator, En);
    Enp1.Add(Enp1_p0n_numerator / denominator, p0n);
    Enp1.Add(Enp1_p1n_numerator / denominator, p1n);

    const double updated_kinetic_energy_density_0 = euler::kineticEnergyDensity(rho_0, p0np1);
    const double updated_kinetic_energy_density_1 = euler::kineticEnergyDensity(rho_1, p1np1);

    const double updated_total_energy_density_0 = updated_kinetic_energy_density_0 + internal_energy_density_0;
    const double updated_total_energy_density_1 = updated_kinetic_energy_density_1 + internal_energy_density_1;

    mfem::Vector updated_euler_dofs_on_cell_0 = euler::constructConservativeState(rho_0, p0np1, updated_total_energy_density_0);
    mfem::Vector updated_euler_dofs_on_cell_1 = euler::constructConservativeState(rho_1, p1np1, updated_total_energy_density_1);

    updated_euler_grid_function_0.SetSubVector(index_of_euler_dofs_on_cell, updated_euler_dofs_on_cell_0);
    updated_euler_grid_function_1.SetSubVector(index_of_euler_dofs_on_cell, updated_euler_dofs_on_cell_1);

    updated_e_field_grid_function.SetSubVector(index_of_e_field_dofs_on_cell, Enp1);
  }

  return std::make_pair(updated_state, updated_field_state);
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
} // namespace mfpic
