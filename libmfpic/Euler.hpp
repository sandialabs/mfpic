#pragma once

#include <mfem/mfem.hpp>

namespace mfpic {

struct Species;

// @brief namespace to encapsulate transformations for the Euler equations, assumes Ideal Gas equation of state
namespace euler {

struct ConservativeVariables {
  enum {
    MASS_DENSITY = 0,
    X_MOMENTUM_DENSITY,
    Y_MOMENTUM_DENSITY,
    Z_MOMENTUM_DENSITY,
    TOTAL_ENERGY_DENSITY,
    NUM_VARS
  };
};

struct PrimitiveVariables {
  enum {
    NUMBER_DENSITY = 0,
    X_BULK_VELOCITY,
    Y_BULK_VELOCITY,
    Z_BULK_VELOCITY,
    TEMPERATURE,
    NUM_VARS
  };
};

double kineticEnergyDensity(const double mass_density, const mfem::Vector& momentum_density);
double pressure(const double number_density, const double temperature);
double pressure(const Species& species, const double internal_energy_density);
double internalEnergyPerUnitMass(const Species& species, const double mass_density, const double pressure);
double internalEnergyDensity(const double internal_energy_density_per_unit_mass, const double mass_density);
double temperature(const double number_density, const double pressure);
double speedOfSound(const Species& species, const double mass_density, const double pressure);

mfem::Vector constructConservativeState(
  const double mass_density,
  const mfem::Vector& momentum_density,
  const double total_energy_density);

mfem::Vector constructPrimitiveState(
  const double number_density,
  const mfem::Vector& bulk_velocity,
  const double temperature);

mfem::Vector convertFromConservativeToPrimitive(const mfem::Vector& conservative_state, const Species& species);
mfem::Vector convertFromPrimitiveToConservative(const mfem::Vector& primitive_state, const Species& species);

double getNumberDensityFromConservativeState(const mfem::Vector& conservative_state, const Species& species);
mfem::Vector getMomentumDensityFromConservativeState(const mfem::Vector& conservative_state);
mfem::Vector getBulkVelocityFromConservativeState(const mfem::Vector& conservative_state);
double getKineticEnergyDensityFromConservativeState(const mfem::Vector& conservative_state);
double getTemperatureFromConservativeState(const mfem::Vector& conservative_state, const Species& species);
double getPressureFromConservativeState(const mfem::Vector& conservative_state, const Species& species);

double getMassDensityFromPrimitiveState(const mfem::Vector& primitive_state, const Species& species);
mfem::Vector getBulkVelocityFromPrimitiveState(const mfem::Vector& primitive_state);
double getPressureFromPrimitiveState(const mfem::Vector& primitive_state);
double getInternalEnergyDensityFromPrimitiveState(const mfem::Vector& primitive_state, const Species& species);

}

}
