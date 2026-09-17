#include <libmfpic/Constants.hpp>
#include <libmfpic/Errors.hpp>
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

double evaluateMaxwellian(const mfem::Vector& primitive_state,
                          const mfem::Vector& velocity,
                          const Species& species)
{
  const double temperature = primitive_state(euler::PrimitiveVariables::TEMPERATURE);
  const double sigma = std::sqrt(constants::boltzmann_constant * temperature / species.mass);

  if ((sigma <= 0.0) || !std::isfinite(sigma)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const int vdim = velocity.Size(); 
  const mfem::Vector bulk_velocity = getBulkVelocityFromPrimitiveState(primitive_state);

  mfem::Vector difference = velocity;
  for (int i=0; i < vdim; ++i )
    difference(i) -= bulk_velocity(i);

  const double inv_sq_sigma = 1.0 / (sigma * sigma);
  const double exponent = inv_sq_sigma * (difference * difference);

  const double norm_base = 1.0 / (std::sqrt(2.0 * M_PI) * sigma);
  const double norm = std::pow(norm_base, static_cast<double>(vdim));

  const double probability_density_function = norm * std::exp(-0.5 * exponent);

  return probability_density_function *
         primitive_state(euler::PrimitiveVariables::NUMBER_DENSITY);
}

double evaluateProductOf1DKappaDistributions(
  const mfem::Vector& primitive_state,
  const mfem::Vector& velocity,
  const double kappa,
  const Species& species)
{
  const double temperature =
    primitive_state(euler::PrimitiveVariables::TEMPERATURE);

  const double number_density =
    primitive_state(euler::PrimitiveVariables::NUMBER_DENSITY);

  if (temperature < 0.0 || species.mass <= 0.0 || kappa <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double v_thermal_squared =
    constants::boltzmann_constant * temperature / species.mass;

  if (!std::isfinite(v_thermal_squared)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double nu = 2.0 * kappa + 1.0;

  const double scale =
    std::sqrt(2.0 * kappa * v_thermal_squared / nu);

  if (scale <= 0.0 || !std::isfinite(scale)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const mfem::Vector bulk_velocity = getBulkVelocityFromPrimitiveState(primitive_state);

  mfem::Vector difference = velocity;
  difference -= bulk_velocity;

  const double scale_squared = scale * scale;
  const double log_norm_1d =
      std::lgamma(0.5 * (nu + 1.0))
    - std::lgamma(0.5 * nu)
    - 0.5 * std::log(nu * M_PI)
    - std::log(scale);

  double log_pdf = std::log(number_density);
  for (int d = 0; d < 3; ++d) {
    const double x = difference(d);
    const double log_shape =
      -0.5 * (nu + 1.0) * std::log1p((x * x) / (nu * scale_squared));

    log_pdf += log_norm_1d + log_shape;
  }
  return std::exp(log_pdf);
}

double evaluateIsotropicKappaDistribution(
    const mfem::Vector& primitive_state,
    const mfem::Vector& velocity,
    const double kappa,
    const Species& species)
{
    const double invalid = std::numeric_limits<double>::quiet_NaN();
    const int velocity_dimensions = velocity.Size(); 

    if (velocity_dimensions != 1 && velocity_dimensions != 3) {
        return invalid;
    }

    const double temperature = primitive_state(PrimitiveVariables::TEMPERATURE);
    const double number_density = primitive_state(PrimitiveVariables::NUMBER_DENSITY);
    const double mass = species.mass;

    if (!std::isfinite(temperature) || temperature <= 0.0 ||
        !std::isfinite(number_density) || number_density < 0.0 ||
        !std::isfinite(mass) || mass <= 0.0 ||
        !std::isfinite(kappa) || kappa <= 1.5) {
        return invalid;
    }

    const mfem::Vector bulk_velocity = getBulkVelocityFromPrimitiveState(primitive_state);
    double thermal_speed = 0.0;
    for (int i = 0; i < velocity_dimensions; ++i) {
        const double fluctuation = velocity(i) - bulk_velocity(i);
        thermal_speed = std::hypot(thermal_speed, fluctuation);
    }

    const double d = static_cast<double>(velocity_dimensions);
    const double half_d = 0.5 * d;
    const double alpha = kappa - 0.5;

    // a^2 = (2*kappa - 3) * k_B*T/m.
    const double log_a_squared =
          std::log(2.0)
        + std::log(kappa - 1.5)
        + std::log(constants::boltzmann_constant)
        + std::log(temperature)
        - std::log(mass);
    const double log_pi = std::log(std::acos(-1.0));
    const double log_normalization = std::lgamma(alpha + half_d) - std::lgamma(alpha) - half_d * (log_pi + log_a_squared);

    // Compute log(1 + |c|^2/a^2) 
    double log_shape_argument = 0.0;
    if (thermal_speed > 0.0) {
        const double log_ratio = 2.0 * std::log(thermal_speed) - log_a_squared;
        log_shape_argument = log_ratio > 0.0 ? log_ratio + std::log1p(std::exp(-log_ratio)) : std::log1p(std::exp(log_ratio));
    }
    const double log_pdf = std::log(number_density) + log_normalization - (alpha + half_d) * log_shape_argument;
    return std::exp(log_pdf);
}

}

}
