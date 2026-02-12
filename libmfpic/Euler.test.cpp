#include <libmfpic/Constants.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/Species.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

TEST(Euler, getTemperatureFromConservativeState) {
  constexpr Species species{.mass = constants::electron_mass};

  constexpr double mass_density = 8.1e-10;
  const mfem::Vector momentum_density{2.5, 3.2, 5.6};
  const double kinetic_energy_density = 0.5 / mass_density * (momentum_density * momentum_density);
  const double total_energy_density = 1.001 * kinetic_energy_density;

  const double number_density = mass_density / species.mass;
  const double internal_energy_density = total_energy_density - kinetic_energy_density;
  const double pressure = internal_energy_density * (species.specific_heat_ratio - 1.);
  const double expected_temperature = pressure / (number_density * constants::boltzmann_constant);

  const mfem::Vector conservative_state = euler::constructConservativeState(mass_density, momentum_density, total_energy_density);

  const double temperature = euler::getTemperatureFromConservativeState(conservative_state, species);
  EXPECT_DOUBLE_EQ(expected_temperature, temperature);
}

TEST(Euler, getPressureFromPrimitiveState) {
  constexpr double number_density = 1.2e16;
  const mfem::Vector bulk_velocity{5.4, 4.7, 8.3};
  constexpr double temperature = 305.1;

  const double expected_pressure = number_density * temperature * constants::boltzmann_constant;

  const mfem::Vector primitive_state = euler::constructPrimitiveState(number_density, bulk_velocity, temperature);
  const double pressure = euler::getPressureFromPrimitiveState(primitive_state);

  EXPECT_DOUBLE_EQ(expected_pressure, pressure);
}

TEST(Euler, getInternalEnergyDensityFromPrimitiveState) {
  constexpr Species species{.mass = constants::electron_mass};

  constexpr double number_density = 1.2e16;
  const mfem::Vector bulk_velocity{5.4, 4.7, 8.3};
  constexpr double temperature = 305.1;

  const double pressure = number_density * temperature * constants::boltzmann_constant;
  const double expected_internal_energy_density = pressure / (species.specific_heat_ratio - 1.);

  const mfem::Vector primitive_state = euler::constructPrimitiveState(number_density, bulk_velocity, temperature);
  const double internal_energy_density = euler::getInternalEnergyDensityFromPrimitiveState(primitive_state, species);

  EXPECT_DOUBLE_EQ(expected_internal_energy_density, internal_energy_density);
}

TEST(Euler, convertFromPrimitiveToConservativeInvertsConvertConservativeToPrimitive) {
  constexpr Species species{.mass = constants::electron_mass};

  constexpr double mass_density = 8.1;
  const mfem::Vector momentum_density{2.5, 3.2, 5.6};
  constexpr double total_energy_density = 9.1;

  const mfem::Vector conservative_state = euler::constructConservativeState(mass_density, momentum_density, total_energy_density);
  const mfem::Vector primitive_state = euler::convertFromConservativeToPrimitive(conservative_state, species);
  const mfem::Vector conservative_state_final = euler::convertFromPrimitiveToConservative(primitive_state, species);

  for (int i = 0; i < conservative_state.Size(); ++i) {
    EXPECT_DOUBLE_EQ(conservative_state[i], conservative_state_final[i]);
  }
}

TEST(Euler, convertFromConservativeToPrimitiveInvertsConvertPrimitiveToConservative) {
  constexpr Species species{.mass = constants::electron_mass};

  constexpr double number_density = 1.2e16;
  const mfem::Vector bulk_velocity{5.4, 4.7, 8.3};
  constexpr double temperature = 305.1;

  const mfem::Vector primitive_state = euler::constructPrimitiveState(number_density, bulk_velocity, temperature);
  const mfem::Vector conservative_state = euler::convertFromPrimitiveToConservative(primitive_state, species);
  const mfem::Vector primitive_state_final = euler::convertFromConservativeToPrimitive(conservative_state, species);

  for (int i = 0; i < primitive_state.Size(); ++i) {
    EXPECT_DOUBLE_EQ(primitive_state[i], primitive_state_final[i]);
  }
}

TEST(Euler, getPressureFromConservativeState) {
  constexpr Species species{.mass = constants::electron_mass};

  constexpr double mass_density = 8.1e-10;
  const mfem::Vector momentum_density{2.5, 3.2, 5.6};
  const double kinetic_energy_density = 0.5 / mass_density * (momentum_density * momentum_density);
  const double total_energy_density = 1.001 * kinetic_energy_density;

  const double internal_energy_density = total_energy_density - kinetic_energy_density;
  const double expected_pressure = internal_energy_density * (species.specific_heat_ratio - 1.);

  const mfem::Vector conservative_state = euler::constructConservativeState(mass_density, momentum_density, total_energy_density);

  const double pressure = euler::getPressureFromConservativeState(conservative_state, species);
  EXPECT_DOUBLE_EQ(expected_pressure, pressure);
}

}
