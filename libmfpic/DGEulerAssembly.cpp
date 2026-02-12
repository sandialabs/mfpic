#include <libmfpic/DGEulerAssembly.hpp>
#include <libmfpic/DGEulerKEIntegrator.hpp>
#include <libmfpic/DGEulerMaxwellSourceIntegrator.hpp>
#include <libmfpic/ElectromagneticFieldsEvaluator.hpp>
#include <libmfpic/Errors.hpp>
#include <libmfpic/LowFidelityState.hpp>
#include <libmfpic/Species.hpp>

#include <mfem/fem/fespace.hpp>
#include <mfem/fem/gridfunc.hpp>
#include <mfem/fem/hyperbolic.hpp>
#include <mfem/linalg/densemat.hpp>

namespace mfpic {

  DGEulerAssembly::DGEulerAssembly(
    mfem::FiniteElementSpace &finite_element_space,
    const Species& species) :
      DGEulerAssembly(CreateDGEulerAssembly_(finite_element_space, species)) {}

  void DGEulerAssembly::computeSources(
    const LowFidelitySpeciesState& species_state,
    const ElectromagneticFieldsEvaluator &field_evaluator,
    mfem::Vector &rhs) const {

    constexpr bool include_energy_source = false;

    const mfem::GridFunction& state_evaluator = species_state.getGridFunction();
    const double charge_over_mass = species_state.getSpecies().charge_over_mass;

    mfem::LinearForm source_form(&getFiniteElementSpace());
    source_form.AddDomainIntegrator(
      new EulerMaxwellSourceIntegrator(state_evaluator, field_evaluator, charge_over_mass, include_energy_source));
    source_form.Assemble();

    rhs += source_form;
  }

  void DGEulerAssembly::computeIntegratedKineticEnergy(const LowFidelitySpeciesState& species_state, mfem::Vector& rhs) const {
    const mfem::GridFunction& state_evaluator = species_state.getGridFunction();

    mfem::LinearForm form(&getFiniteElementSpace());
    form.AddDomainIntegrator(new EulerKineticEnergyIntegrator(state_evaluator));
    form.Assemble();

    rhs += form;
  }

} // namespace mfpic
