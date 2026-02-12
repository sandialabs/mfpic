#include <libmfpic/Constants.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/Species.hpp>

namespace mfpic {

namespace euler {

double kineticEnergyDensity(const double mass_density, const mfem::Vector& momentum_density) {
  return 0.5 / mass_density * (momentum_density * momentum_density);
}

double pressure(const double number_density, const double temperature) {
  return number_density * temperature * constants::boltzmann_constant;
}

double pressure(const Species& species, const double internal_energy_density) {
  return (species.specific_heat_ratio - 1.) * internal_energy_density;
}

double internalEnergyPerUnitMass(const Species& species, const double mass_density, const double pressure) {
  return pressure / ((species.specific_heat_ratio - 1.) * mass_density);
}

double internalEnergyDensity(const double internal_energy_density_per_unit_mass, const double mass_density) {
  return internal_energy_density_per_unit_mass * mass_density;
}

double temperature(const double number_density, const double pressure) {
  return pressure / (number_density * constants::boltzmann_constant);
}

double speedOfSound(const Species& species, const double mass_density, const double pressure) {
  return std::sqrt(species.specific_heat_ratio * pressure / mass_density);
}

mfem::Vector constructConservativeState(
  const double mass_density,
  const mfem::Vector& momentum_density,
  const double total_energy_density)
{
  mfem::Vector conservative_state(ConservativeVariables::NUM_VARS);
  conservative_state[ConservativeVariables::MASS_DENSITY] = mass_density;
  conservative_state[ConservativeVariables::X_MOMENTUM_DENSITY] = momentum_density[0];
  conservative_state[ConservativeVariables::Y_MOMENTUM_DENSITY] = momentum_density[1];
  conservative_state[ConservativeVariables::Z_MOMENTUM_DENSITY] = momentum_density[2];
  conservative_state[ConservativeVariables::TOTAL_ENERGY_DENSITY] = total_energy_density;

  return conservative_state;
}

mfem::Vector constructPrimitiveState(
  const double number_density,
  const mfem::Vector& bulk_velocity,
  const double temperature)
{
  mfem::Vector primitive_state(PrimitiveVariables::NUM_VARS);
  primitive_state[PrimitiveVariables::NUMBER_DENSITY] = number_density;
  primitive_state[PrimitiveVariables::X_BULK_VELOCITY] = bulk_velocity[0];
  primitive_state[PrimitiveVariables::Y_BULK_VELOCITY] = bulk_velocity[1];
  primitive_state[PrimitiveVariables::Z_BULK_VELOCITY] = bulk_velocity[2];
  primitive_state[PrimitiveVariables::TEMPERATURE] = temperature;

  return primitive_state;
}

mfem::Vector convertFromConservativeToPrimitive(const mfem::Vector& conservative_state, const Species& species) {
  const double number_density = getNumberDensityFromConservativeState(conservative_state, species);
  const mfem::Vector bulk_velocity = getBulkVelocityFromConservativeState(conservative_state);
  const double temperature = getTemperatureFromConservativeState(conservative_state, species);
  return constructPrimitiveState(number_density, bulk_velocity, temperature);
}

mfem::Vector convertFromPrimitiveToConservative(const mfem::Vector& primitive_state, const Species& species) {
  const double mass_density = getMassDensityFromPrimitiveState(primitive_state, species);

  mfem::Vector momentum_density = getBulkVelocityFromPrimitiveState(primitive_state);
  momentum_density *= mass_density;

  const double kinetic_energy_density = kineticEnergyDensity(mass_density, momentum_density);
  const double internal_energy_density = getInternalEnergyDensityFromPrimitiveState(primitive_state, species);
  const double total_energy_density = kinetic_energy_density + internal_energy_density;

  return constructConservativeState(mass_density, momentum_density, total_energy_density);
}

double getNumberDensityFromConservativeState(const mfem::Vector& conservative_state, const Species& species) {
  return conservative_state[ConservativeVariables::MASS_DENSITY] / species.mass;
}

mfem::Vector getMomentumDensityFromConservativeState(const mfem::Vector& conservative_state) {
  mfem::Vector momentum_density{
    conservative_state[ConservativeVariables::X_MOMENTUM_DENSITY],
    conservative_state[ConservativeVariables::Y_MOMENTUM_DENSITY],
    conservative_state[ConservativeVariables::Z_MOMENTUM_DENSITY]};

  return momentum_density;
}

mfem::Vector getBulkVelocityFromConservativeState(const mfem::Vector& conservative_state) {
  const double mass_density = conservative_state[ConservativeVariables::MASS_DENSITY];
  const mfem::Vector bulk_velocity{
    conservative_state[ConservativeVariables::X_MOMENTUM_DENSITY] / mass_density,
    conservative_state[ConservativeVariables::Y_MOMENTUM_DENSITY] / mass_density,
    conservative_state[ConservativeVariables::Z_MOMENTUM_DENSITY] / mass_density};

  return bulk_velocity;
}

double getKineticEnergyDensityFromConservativeState(const mfem::Vector& conservative_state) {
  const double mass_density = conservative_state[ConservativeVariables::MASS_DENSITY];
  const mfem::Vector momentum_density = getMomentumDensityFromConservativeState(conservative_state);
  return kineticEnergyDensity(mass_density, momentum_density);
}

double getTemperatureFromConservativeState(const mfem::Vector& conservative_state, const Species& species) {
  const double number_density = getNumberDensityFromConservativeState(conservative_state, species);
  const double pressure_value = getPressureFromConservativeState(conservative_state, species);
  const double temperature_value = temperature(number_density, pressure_value);
  return temperature_value;
}

double getPressureFromConservativeState(const mfem::Vector& conservative_state, const Species& species) {
  const double total_energy_density = conservative_state[ConservativeVariables::TOTAL_ENERGY_DENSITY];
  const double kinetic_energy_density = getKineticEnergyDensityFromConservativeState(conservative_state);
  const double internal_energy_density = total_energy_density - kinetic_energy_density;
  const double pressure_value = pressure(species, internal_energy_density);
  return pressure_value;
}

double getMassDensityFromPrimitiveState(const mfem::Vector& primitive_state, const Species& species) {
  return primitive_state[PrimitiveVariables::NUMBER_DENSITY] * species.mass;
}

mfem::Vector getBulkVelocityFromPrimitiveState(const mfem::Vector& primitive_state) {
  const mfem::Vector bulk_velocity {
    primitive_state[PrimitiveVariables::X_BULK_VELOCITY],
    primitive_state[PrimitiveVariables::Y_BULK_VELOCITY],
    primitive_state[PrimitiveVariables::Z_BULK_VELOCITY]};
  return bulk_velocity;
}

double getPressureFromPrimitiveState(const mfem::Vector& primitive_state) {
  const double number_density = primitive_state[PrimitiveVariables::NUMBER_DENSITY];
  const double temperature = primitive_state[PrimitiveVariables::TEMPERATURE];
  return pressure(number_density, temperature);
}

double getInternalEnergyDensityFromPrimitiveState(const mfem::Vector& primitive_state, const Species& species) {
  const double mass_density = getMassDensityFromPrimitiveState(primitive_state, species);
  const double pressure = getPressureFromPrimitiveState(primitive_state);
  const double internal_energy_density_per_unit_mass = internalEnergyPerUnitMass(species, mass_density, pressure);
  const double internal_energy_density = internalEnergyDensity(internal_energy_density_per_unit_mass, mass_density);
  return internal_energy_density;
}

}

}
