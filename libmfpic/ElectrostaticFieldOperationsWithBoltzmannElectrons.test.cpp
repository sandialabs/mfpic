#include <libmfpic/Constants.hpp>
#include <libmfpic/DirichletBoundaryConditionsConstant.hpp>
#include <libmfpic/ElectrostaticFieldOperationsWithBoltzmannElectrons.hpp>
#include <libmfpic/ElectrostaticFieldState.hpp>
#include <libmfpic/MeshUtilities.hpp>
#include <libmfpic/Pinning.hpp>

#include <gtest/gtest.h>
#include <mfem/mfem.hpp>

using namespace mfpic;

TEST(ElectrostaticFieldOperationsWithBoltzmannElectrons, ZeroChargeAndZeroReferenceDensityWithNaturalBoundariesGivesZeroPotential) {
  constexpr int num_elems = 10;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);
  constexpr int hgrad_order = 1;
  Discretization es_discretization(&mesh, hgrad_order);
  constexpr double reference_number_density = 0.0;
  constexpr double temperature = 1.0;
  mfpic::ElectrostaticFieldOperationsWithBoltzmannElectrons es_field_operations(
    es_discretization,
    std::make_unique<Pinning>(),
    reference_number_density,
    temperature
  );

  ElectrostaticFieldState es_field_state(es_discretization);
  IntegratedCharge integrated_charge(es_discretization);
  integrated_charge.setIntegratedChargeValue(0.);
  es_field_operations.fieldSolve(es_field_state, integrated_charge);

  mfem::GridFunction potential = es_field_state.getPotential();
  mfem::ConstantCoefficient exact_potential(0.);
  const double l2_error = potential.ComputeL2Error(exact_potential);
  constexpr double expected_error = 0.;
  EXPECT_EQ(expected_error, l2_error);
}
