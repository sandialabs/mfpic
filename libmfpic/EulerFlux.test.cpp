#include <libmfpic/Constants.hpp>
#include <libmfpic/EulerFlux.hpp>
#include <libmfpic/Species.hpp>

#include <gtest/gtest.h>
#include <mfem/fem/eltrans.hpp>
#include <mfem/linalg/densemat.hpp>

namespace {

using namespace mfpic;

auto run_flux_test = [] (const int spatial_dim) {

  constexpr double specific_heat_ratio = 1.4;
  constexpr Species species{.mass = constants::electron_mass, .specific_heat_ratio = specific_heat_ratio};

  constexpr double mass_density = 8.1e-10;
  const double number_density = mass_density / species.mass;
  const mfem::Vector momentum_density{2.5, 3.2, 5.6};
  const mfem::Vector velocity {momentum_density[0] / mass_density,
                               momentum_density[1] / mass_density,
                               momentum_density[2] / mass_density};
  constexpr double temperature = 505.0;
  const double pressure = number_density * temperature * constants::boltzmann_constant; 
  const double kinetic_energy_density = 0.5 / mass_density * (momentum_density * momentum_density);
  const double internal_energy_density = pressure / (species.specific_heat_ratio - 1.);
  const double total_energy_density = kinetic_energy_density + internal_energy_density;

  const mfem::Vector conservative_state = euler::constructConservativeState(mass_density, momentum_density, total_energy_density);

  mfem::IsoparametricTransformation tr; // does nothing
  mfem::DenseMatrix flux(5,spatial_dim), expected_flux(5,spatial_dim);

  EulerFlux flux_function(spatial_dim,species);

  flux_function.ComputeFlux(conservative_state, tr, flux);

  switch (spatial_dim) { 
  case 3:
    expected_flux(0,2) = momentum_density(2);
    expected_flux(1,2) = momentum_density(0) * velocity(2);
    expected_flux(2,2) = momentum_density(1) * velocity(2);
    expected_flux(3,2) = momentum_density(2) * velocity(2) + pressure;
    expected_flux(4,2) = (total_energy_density + pressure) * velocity(2);
    [[fallthrough]];
  case 2:
    expected_flux(0,1) = momentum_density(1);
    expected_flux(1,1) = momentum_density(0) * velocity(1);
    expected_flux(2,1) = momentum_density(1) * velocity(1) + pressure;
    expected_flux(3,1) = momentum_density(2) * velocity(1);
    expected_flux(4,1) = (total_energy_density + pressure) * velocity(1);
    [[fallthrough]];
  case 1:
    expected_flux(0,0) = momentum_density(0);
    expected_flux(1,0) = momentum_density(0) * velocity(0) + pressure;
    expected_flux(2,0) = momentum_density(1) * velocity(0);
    expected_flux(3,0) = momentum_density(2) * velocity(0);
    expected_flux(4,0) = (total_energy_density + pressure) * velocity(0);
    break;
  }

  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < spatial_dim; ++j)
      EXPECT_DOUBLE_EQ(expected_flux(i,j), flux(i,j)) << "Component: " << i << " Dim: " << j;
};

auto run_flux_dot_n_test = [] (const int spatial_dim, const mfem::Vector normal) {

  ASSERT_TRUE(normal.Size() == spatial_dim) << "To run this test, normal must match spatial dimension!";

  constexpr double specific_heat_ratio = 1.4;
  constexpr Species species{.mass = constants::electron_mass, .specific_heat_ratio = specific_heat_ratio};

  constexpr double mass_density = 8.1e-10;
  const double number_density = mass_density / species.mass;
  const mfem::Vector momentum_density{2.5, 3.2, 5.6};
  const mfem::Vector velocity {momentum_density[0] / mass_density,
                               momentum_density[1] / mass_density,
                               momentum_density[2] / mass_density};
  constexpr double temperature = 505.0;
  const double pressure = number_density * temperature * constants::boltzmann_constant; 
  const double kinetic_energy_density = 0.5 / mass_density * (momentum_density * momentum_density);
  const double internal_energy_density = pressure / (species.specific_heat_ratio - 1.);
  const double total_energy_density = kinetic_energy_density + internal_energy_density;

  const mfem::Vector conservative_state = euler::constructConservativeState(mass_density, momentum_density, total_energy_density);

  mfem::FaceElementTransformations tr; // does nothing
  mfem::Vector flux_dot_n(5), expected_flux_dot_n(5);
  mfem::DenseMatrix expected_flux(5,spatial_dim);

  EulerFlux flux_function(spatial_dim,species);

  flux_function.ComputeFluxDotN(conservative_state, normal, tr, flux_dot_n);

  switch (spatial_dim) { 
  case 3:
    expected_flux(0,2) = momentum_density(2);
    expected_flux(1,2) = momentum_density(0) * velocity(2);
    expected_flux(2,2) = momentum_density(1) * velocity(2);
    expected_flux(3,2) = momentum_density(2) * velocity(2) + pressure;
    expected_flux(4,2) = (total_energy_density + pressure) * velocity(2);
    [[fallthrough]];
  case 2:
    expected_flux(0,1) = momentum_density(1);
    expected_flux(1,1) = momentum_density(0) * velocity(1);
    expected_flux(2,1) = momentum_density(1) * velocity(1) + pressure;
    expected_flux(3,1) = momentum_density(2) * velocity(1);
    expected_flux(4,1) = (total_energy_density + pressure) * velocity(1);
    [[fallthrough]];
  case 1:
    expected_flux(0,0) = momentum_density(0);
    expected_flux(1,0) = momentum_density(0) * velocity(0) + pressure;
    expected_flux(2,0) = momentum_density(1) * velocity(0);
    expected_flux(3,0) = momentum_density(2) * velocity(0);
    expected_flux(4,0) = (total_energy_density + pressure) * velocity(0);
    break;
  }

  expected_flux.Mult(normal,expected_flux_dot_n);

  for (int i = 0; i < 5; ++i)
    EXPECT_DOUBLE_EQ(expected_flux_dot_n(i), flux_dot_n(i)) << "Component: " << i;
};

TEST(EulerFlux, ComputeFlux3D) {
  run_flux_test(3);
}
TEST(EulerFlux, ComputeFlux2D) {
  run_flux_test(2);
}
TEST(EulerFlux, ComputeFlux1D) {
  run_flux_test(1);
}

TEST(EulerFlux, ComputeFluxDotN3D) {
  mfem::Vector x{1.,0.,0.};
  mfem::Vector y{0.,1.,0.};
  mfem::Vector z{0.,0.,1.};
  mfem::Vector all{.2345,-.123,.789}; // does not need to be normalized (and usually isn't internally)
  run_flux_dot_n_test(3,x);
  run_flux_dot_n_test(3,y);
  run_flux_dot_n_test(3,z);
  run_flux_dot_n_test(3,all);
}
TEST(EulerFlux, ComputeFluxDotN2D) {
  mfem::Vector x{1.,0.};
  mfem::Vector y{0.,1.};
  mfem::Vector all{.2345,-.123}; // does not need to be normalized (and usually isn't internally)
  run_flux_dot_n_test(2,x);
  run_flux_dot_n_test(2,y);
  run_flux_dot_n_test(2,all);
}
TEST(EulerFlux, ComputeFluxDotN1D) {
  mfem::Vector x{1.};
  run_flux_dot_n_test(1,x);
}
}
