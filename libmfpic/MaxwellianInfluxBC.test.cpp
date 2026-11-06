#include <libmfpic/Constants.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/KineticFluxBC.hpp>
#include <libmfpic/MaxwellianInfluxBC.hpp>

#include <gtest/gtest.h>
#include <mfem/mfem.hpp>

namespace {

using namespace mfpic;

const Species electron_species{
  .charge = -constants::elementary_charge,
  .mass = constants::electron_mass,
  .specific_heat_ratio = 5./3.
};

mfem::Vector evalMaxwellianInfluxFlux(
  int spatial_dim,
  double number_flux_density,
  double temperature,
  const mfem::Vector& conservative_state,
  const mfem::Vector& normal
) {
  constexpr int num_equations = 5;
  DummyFluxFunction dummy_flux(num_equations, spatial_dim);
  MaxwellianInfluxNumericalFlux flux(dummy_flux, electron_species, number_flux_density, temperature);
  mfem::FaceElementTransformations dummy_transform;
  mfem::Vector f_dot_n(num_equations);
  flux.Eval(conservative_state, conservative_state, normal, dummy_transform, f_dot_n);
  return f_dot_n;
}

double reservoirNumberDensity(double number_flux_density, double temperature) {
  return number_flux_density
    / std::sqrt(constants::boltzmann_constant * temperature / (2. * M_PI * electron_species.mass));
}

TEST(MaxwellianInfluxBC, InfluxIntoEmptyDomainGivesHalfRangeMoments) {
  constexpr double number_flux_density = 2.3e21;
  constexpr double temperature = 11600.0;
  constexpr double kT = constants::boltzmann_constant * temperature;
  const double half_range_pressure = reservoirNumberDensity(number_flux_density, temperature) * kT / 2.;
  constexpr double vacuum_state_number_density = 1.0e-30;
  const mfem::Vector vacuum_state = euler::convertFromPrimitiveToConservative(
    euler::constructPrimitiveState(vacuum_state_number_density, mfem::Vector{0., 0., 0.}, temperature), electron_species);

  for (int spatial_dim = 1; spatial_dim <= 3; spatial_dim++) {
    mfem::Vector normal{.65, -.2, .11};
    for (int idim = spatial_dim; idim < 3; idim++) {
      normal(idim) = 0.0;
    }
    mfem::Vector unit_normal = normal;
    unit_normal /= normal.Norml2();

    mfem::Vector f_dot_n_expected(5);
    f_dot_n_expected = 0.;
    f_dot_n_expected(0) = -electron_species.mass * number_flux_density;
    for (int d = 0; d < spatial_dim; ++d)
      f_dot_n_expected(1 + d) = half_range_pressure * unit_normal(d);
    f_dot_n_expected(4) = -2. * kT * number_flux_density;
    f_dot_n_expected *= normal.Norml2();

    const mfem::Vector f_dot_n = evalMaxwellianInfluxFlux(spatial_dim, number_flux_density, temperature, vacuum_state, normal);
    for (int i = 0; i < 5; ++i) {
      EXPECT_DOUBLE_EQ(f_dot_n(i), f_dot_n_expected(i)) << spatial_dim << " " << i;
    }
  }
}

TEST(MaxwellianInfluxBC, InfluxInEquilibriumWithReservoirHasZeroNetFlux) {
  constexpr int spatial_dim = 3;
  constexpr double number_flux_density = 2.3e21;
  constexpr double temperature = 11600.0;
  const mfem::Vector reservoir_state = euler::convertFromPrimitiveToConservative(
    euler::constructPrimitiveState(reservoirNumberDensity(number_flux_density, temperature), mfem::Vector{0., 0., 0.}, temperature),
    electron_species);
  const double pressure = euler::getPressureFromConservativeState(reservoir_state, electron_species);

  const mfem::Vector normal {.65, -.2, .11};
  const mfem::Vector f_dot_n = evalMaxwellianInfluxFlux(spatial_dim, number_flux_density, temperature, reservoir_state, normal);

  constexpr double relative_tolerance_for_zero_comparisons = 1e-12;
  mfem::Vector unit_normal = normal;
  unit_normal /= normal.Norml2();
  const double one_way_mass_flux = electron_species.mass * number_flux_density * normal.Norml2();
  const double one_way_energy_flux = 2. * constants::boltzmann_constant * temperature * number_flux_density * normal.Norml2();
  EXPECT_NEAR(f_dot_n(0), 0., relative_tolerance_for_zero_comparisons * one_way_mass_flux);
  for (int d = 0; d < spatial_dim; ++d) {
    EXPECT_DOUBLE_EQ(f_dot_n(1 + d), pressure * normal(d)) << d;
  }
  EXPECT_NEAR(f_dot_n(4), 0., relative_tolerance_for_zero_comparisons * one_way_energy_flux);
}

TEST(MaxwellianInfluxBC, ZeroInfluxMatchesKineticFlux) {
  constexpr int spatial_dim = 2;
  const mfem::Vector conservative_state = euler::convertFromPrimitiveToConservative(
    euler::constructPrimitiveState(1.567e16, mfem::Vector{1.e5, -3.e4, 2.e4}, 293.), electron_species);
  const mfem::Vector normal {.3, -1.2, 0.};

  constexpr int num_equations = 5;
  DummyFluxFunction dummy_flux(num_equations, spatial_dim);
  KineticFluxNumericalFlux absorbing_flux(dummy_flux, electron_species);
  mfem::FaceElementTransformations dummy_transform;
  mfem::Vector f_dot_n_absorbing(num_equations);
  absorbing_flux.Eval(conservative_state, conservative_state, normal, dummy_transform, f_dot_n_absorbing);

  constexpr double number_flux_density = 1.0e-30;
  constexpr double temperature = 300.0;
  const mfem::Vector f_dot_n = evalMaxwellianInfluxFlux(spatial_dim, number_flux_density, temperature, conservative_state, normal);
  for (int i = 0; i < 5; ++i) {
    EXPECT_DOUBLE_EQ(f_dot_n(i), f_dot_n_absorbing(i)) << i;
  }
}

} // namespace
