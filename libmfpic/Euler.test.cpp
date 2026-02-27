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

TEST(Maxwellian, MaxwellianValueAtMean)
{
  constexpr double number_density = 1.2e16;
  const mfem::Vector bulk_velocity{5.4, 4.7, 8.3};
  constexpr double temperature = 305.1;
  constexpr Species species{.mass = constants::electron_mass};
  mfem::Vector prim(5);
  prim(euler::PrimitiveVariables::NUMBER_DENSITY) = number_density;
  prim(euler::PrimitiveVariables::X_BULK_VELOCITY) = bulk_velocity(0);
  prim(euler::PrimitiveVariables::Y_BULK_VELOCITY) = bulk_velocity(1);
  prim(euler::PrimitiveVariables::Z_BULK_VELOCITY) = bulk_velocity(2);
  prim(euler::PrimitiveVariables::TEMPERATURE) = temperature;

  const double sigma = std::sqrt(constants::boltzmann_constant * temperature / species.mass);

  for (int i = 0; i < 3; ++i)
  {
    double val_at_mean = euler::evaluateMaxwellian(prim, bulk_velocity, species,i+1);
    double expected_at_mean =
        number_density / std::pow(std::sqrt(2.0 * M_PI) * sigma, i+1);
    EXPECT_NEAR(val_at_mean, expected_at_mean, expected_at_mean * 1e-12);
  }
}

TEST(Maxwellian, MaxwellianShiftInV)
{
  constexpr double number_density = 1.2e16;
  const mfem::Vector bulk_velocity{5.4, 4.7, 8.3};
  constexpr double temperature = 305.1;
  constexpr Species species{.mass = constants::electron_mass};
  mfem::Vector prim(5);
  prim(euler::PrimitiveVariables::NUMBER_DENSITY) = number_density;
  prim(euler::PrimitiveVariables::X_BULK_VELOCITY) = bulk_velocity(0);
  prim(euler::PrimitiveVariables::Y_BULK_VELOCITY) = bulk_velocity(1);
  prim(euler::PrimitiveVariables::Z_BULK_VELOCITY) = bulk_velocity(2);
  prim(euler::PrimitiveVariables::TEMPERATURE) = temperature;

  const double sigma = std::sqrt(constants::boltzmann_constant * temperature / species.mass);

  const double val_at_mean = euler::evaluateMaxwellian(prim, bulk_velocity, species,3);

  const double a = 0.7 * sigma;
  mfem::Vector v_shift(3);
  v_shift(0) = bulk_velocity(0) + a; v_shift(1) = bulk_velocity(1); v_shift(2) = bulk_velocity(2);

  const double val_shift = euler::evaluateMaxwellian(prim, v_shift, species,3);
  const double expected_ratio = std::exp(-0.5 * (a*a) / (sigma*sigma));
  EXPECT_NEAR(val_shift / val_at_mean, expected_ratio, 1e-12);
}

TEST(Maxwellian, MaxwellianIntegratesToOnein1D)
{
  constexpr double number_density = 1.0;
  const mfem::Vector bulk_velocity{5.4, 4.7, 8.3};
  constexpr double temperature = 305.1;
  constexpr Species species{.mass = constants::electron_mass};
  mfem::Vector prim(5);
  prim(euler::PrimitiveVariables::NUMBER_DENSITY) = number_density;
  prim(euler::PrimitiveVariables::X_BULK_VELOCITY) = bulk_velocity(0);
  prim(euler::PrimitiveVariables::Y_BULK_VELOCITY) = bulk_velocity(1);
  prim(euler::PrimitiveVariables::Z_BULK_VELOCITY) = bulk_velocity(2);
  prim(euler::PrimitiveVariables::TEMPERATURE) = temperature;

  const double sigma = std::sqrt(constants::boltzmann_constant * temperature / species.mass);
  const double L = 6.0 * sigma;
  const int N = 41; 
  const double dv = 2.0 * L / (N - 1);
  double integral = 0.0;
  mfem::Vector v(3);
  for (int i = 0; i < N; ++i) {
    v(0) = bulk_velocity(0) + (-L + i * dv);
    v(1) = 0;
    v(2) = 0;
    integral += euler::evaluateMaxwellian(prim, v, species,1);
  }
  integral *= (dv);
  EXPECT_NEAR(integral, 1.0, 1e-5);
}

TEST(Maxwellian, MaxwellianIntegratesToOnein2D)
{
  constexpr double number_density = 1.0;
  const mfem::Vector bulk_velocity{5.4, 4.7, 8.3};
  constexpr double temperature = 305.1;
  constexpr Species species{.mass = constants::electron_mass};
  mfem::Vector prim(5);
  prim(euler::PrimitiveVariables::NUMBER_DENSITY) = number_density;
  prim(euler::PrimitiveVariables::X_BULK_VELOCITY) = bulk_velocity(0);
  prim(euler::PrimitiveVariables::Y_BULK_VELOCITY) = bulk_velocity(1);
  prim(euler::PrimitiveVariables::Z_BULK_VELOCITY) = bulk_velocity(2);
  prim(euler::PrimitiveVariables::TEMPERATURE) = temperature;

  const double sigma = std::sqrt(constants::boltzmann_constant * temperature / species.mass);
  const double L = 6.0 * sigma;
  const int N = 41; 
  const double dv = 2.0 * L / (N - 1);
  double integral = 0.0;
  mfem::Vector v(3);
  for (int i = 0; i < N; ++i) {
    v(0) = bulk_velocity(0) + (-L + i * dv);
    for (int j = 0; j < N; ++j) {
      v(1) = bulk_velocity(1) + (-L + j * dv);
      v(2) = 0;
      integral += euler::evaluateMaxwellian(prim, v, species,2);
    }
  }
  integral *= (dv * dv);
  EXPECT_NEAR(integral, 1.0, 1e-5);
}

TEST(Maxwellian, MaxwellianIntegratesToOnein3D)
{
  constexpr double number_density = 1.0;
  const mfem::Vector bulk_velocity{5.4, 4.7, 8.3};
  constexpr double temperature = 305.1;
  constexpr Species species{.mass = constants::electron_mass};
  mfem::Vector prim(5);
  prim(euler::PrimitiveVariables::NUMBER_DENSITY) = number_density;
  prim(euler::PrimitiveVariables::X_BULK_VELOCITY) = bulk_velocity(0);
  prim(euler::PrimitiveVariables::Y_BULK_VELOCITY) = bulk_velocity(1);
  prim(euler::PrimitiveVariables::Z_BULK_VELOCITY) = bulk_velocity(2);
  prim(euler::PrimitiveVariables::TEMPERATURE) = temperature;

  const double sigma = std::sqrt(constants::boltzmann_constant * temperature / species.mass);
  const double L = 6.0 * sigma;
  const int N = 41; 
  const double dv = 2.0 * L / (N - 1);
  double integral = 0.0;
  mfem::Vector v(3);
  for (int i = 0; i < N; ++i) {
    v(0) = bulk_velocity(0) + (-L + i * dv);
    for (int j = 0; j < N; ++j) {
      v(1) = bulk_velocity(1) + (-L + j * dv);
      for (int k = 0; k < N; ++k) {
        v(2) = bulk_velocity(2) + (-L + k * dv);
        integral += euler::evaluateMaxwellian(prim, v, species,3);
      }
    }
  }
  integral *= (dv * dv * dv);
  EXPECT_NEAR(integral, 1.0, 1e-5);
}

}
