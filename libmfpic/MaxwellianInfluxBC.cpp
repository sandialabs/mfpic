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
  const double outflow_max_speed = KineticFluxNumericalFlux::Eval(conservative_state, not_used, normal, transformation, flux_dot_n);

  mfem::Vector inward_normal = normal;
  inward_normal *= -1.0;

  const double mean_speed = std::sqrt(8.0 * constants::boltzmann_constant * temperature_ / (M_PI * species.mass));
  const double reservoir_number_density = 4.0 * number_flux_density_ / mean_speed;
  const mfem::Vector reservoir_bulk_velocity{0.0, 0.0, 0.0};
  const mfem::Vector reservoir_primitive_state = constructPrimitiveState(reservoir_number_density, reservoir_bulk_velocity, temperature_);
  mfem::Vector influx_dot_n(ConservativeVariables::NUM_VARS);
  const double inflow_max_speed = KineticFluxNumericalFlux::Eval(
    convertFromPrimitiveToConservative(reservoir_primitive_state, species),
    not_used,
    inward_normal,
    transformation,
    influx_dot_n
  );

  flux_dot_n -= influx_dot_n;

  return std::max(outflow_max_speed, inflow_max_speed);
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
