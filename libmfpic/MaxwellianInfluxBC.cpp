#include <libmfpic/Constants.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/MaxwellianInfluxBC.hpp>

namespace mfpic {

MaxwellianInfluxNumericalFlux::MaxwellianInfluxNumericalFlux(
  const mfem::FluxFunction& flux_function,
  const Species species,
  double number_flux_density,
  double temperature
) :
  KineticFluxNumericalFlux(flux_function, species),
  number_flux_density_(number_flux_density),
  temperature_(temperature)
{}

double MaxwellianInfluxNumericalFlux::Eval(const mfem::Vector &conservative_state,
                                            const mfem::Vector &not_used,
                                            const mfem::Vector &normal,
                                            mfem::FaceElementTransformations &transformation,
                                            mfem::Vector &flux_dot_n) const
{
  using namespace euler;
  const double max_speed = KineticFluxNumericalFlux::Eval(conservative_state, not_used, normal, transformation, flux_dot_n);

  const double d_area = normal.Norml2();
  flux_dot_n(ConservativeVariables::MASS_DENSITY) -= d_area * species.mass * number_flux_density_;

  const double mean_speed = std::sqrt(8.0 * constants::boltzmann_constant * temperature_ / (M_PI * species.mass));
  const double reservoir_number_density = 4.0 * number_flux_density_ / mean_speed;
  const double pressure_from_incoming_particles = reservoir_number_density * constants::boltzmann_constant * temperature_ / 2.0;
  switch (fluxFunction.dim) {
  case 3:
    flux_dot_n(ConservativeVariables::Z_MOMENTUM_DENSITY) += pressure_from_incoming_particles * normal(2);
    [[fallthrough]];
  case 2:
    flux_dot_n(ConservativeVariables::Y_MOMENTUM_DENSITY) += pressure_from_incoming_particles * normal(1);
    [[fallthrough]];
  case 1:
    flux_dot_n(ConservativeVariables::X_MOMENTUM_DENSITY) += pressure_from_incoming_particles * normal(0);
    break;
  }

  const double mean_energy_per_incoming_particle = 2.0 * constants::boltzmann_constant * temperature_;
  flux_dot_n(ConservativeVariables::TOTAL_ENERGY_DENSITY) -= d_area * mean_energy_per_incoming_particle * number_flux_density_;

  return max_speed;
}

MaxwellianInfluxBC::MaxwellianInfluxBC(
  int boundary_attribute,
  const mfem::Mesh& mesh,
  const Species& species,
  double number_flux_density_in,
  double temperature_in
) :
  KineticFluxBC(boundary_attribute, mesh, species),
  number_flux_density(number_flux_density_in),
  temperature(temperature_in)
{};

std::unique_ptr<KineticFluxNumericalFlux> MaxwellianInfluxBC::makeNumericalFlux_(
  const mfem::FluxFunction & flux_function
) const {
  return std::make_unique<MaxwellianInfluxNumericalFlux>(flux_function, species, number_flux_density, temperature);
}

} // namespace mfpic
