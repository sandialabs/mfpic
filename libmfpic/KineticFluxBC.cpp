#include <libmfpic/Constants.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/KineticFluxBC.hpp>

namespace mfpic {

  double KineticFluxFluxFunction::ComputeFluxDotN(const mfem::Vector &conservative_state, 
                                                  const mfem::Vector &normal,
                                                  mfem::FaceElementTransformations &,
                                                  mfem::Vector &flux_dot_n) const
  {
    using namespace euler;
    const double pressure = getPressureFromConservativeState(conservative_state, species_);
    const double temperature = getTemperatureFromConservativeState(conservative_state, species_);
    const mfem::Vector momentum = getMomentumDensityFromConservativeState(conservative_state);
    const mfem::Vector velocity = getBulkVelocityFromConservativeState(conservative_state);

    const double d_area = normal.Norml2(); 
    mfem::Vector unit_normal = normal;
    unit_normal /= d_area;

    // Check whether the solution is physical only in debug mode
    MFEM_ASSERT(conservative_state[ConservativeVariables::MASS_DENSITY] >= 0, "Negative Density");
    MFEM_ASSERT(pressure >= 0, "Negative Pressure");
    MFEM_ASSERT(temperature >= 0, "Negative Temperature");
    MFEM_ASSERT(conservative_state[ConservativeVariables::TOTAL_ENERGY_DENSITY] >= 0, "Negative Total Energy Density");

    double normal_velocity(0.);
    for (int d = 0; d < dim; d++)
      normal_velocity += velocity(d) * unit_normal(d);
    const double most_probable_speed = std::sqrt(2. * constants::boltzmann_constant * temperature / species_.mass);
    const double reduced_velocity = normal_velocity / most_probable_speed;
    const double sqrt_pi = std::sqrt(M_PI);
    const double exp_term = std::exp(-reduced_velocity * reduced_velocity);
    const double erf_term = std::erf(reduced_velocity);

    auto func_erf = [&] (const double & sign) {
      return 1. / 2. * (1. + sign * erf_term);
    };
    auto func_exp = [&] (const double & sign) {
      return sign * most_probable_speed / (2. * sqrt_pi) * exp_term;
    };
    const double plus = 1.; // assume only particle absorption
    const double a_plus = func_erf(plus);
    const double b_plus = func_exp(plus);

    flux_dot_n(ConservativeVariables::MASS_DENSITY) = 
      conservative_state[ConservativeVariables::MASS_DENSITY] *
      ( normal_velocity * a_plus + b_plus );

    flux_dot_n(ConservativeVariables::X_MOMENTUM_DENSITY) = 
      flux_dot_n(ConservativeVariables::MASS_DENSITY) * velocity(0);

    flux_dot_n(ConservativeVariables::Y_MOMENTUM_DENSITY) = 
      flux_dot_n(ConservativeVariables::MASS_DENSITY) * velocity(1);

    flux_dot_n(ConservativeVariables::Z_MOMENTUM_DENSITY) = 
      flux_dot_n(ConservativeVariables::MASS_DENSITY) * velocity(2);

    switch (dim) {
    case 3:
      flux_dot_n(ConservativeVariables::Z_MOMENTUM_DENSITY) += pressure * unit_normal(2) * a_plus;
      [[fallthrough]];
    case 2:
      flux_dot_n(ConservativeVariables::Y_MOMENTUM_DENSITY) += pressure * unit_normal(1) * a_plus;
      [[fallthrough]];
    case 1:
      flux_dot_n(ConservativeVariables::X_MOMENTUM_DENSITY) += pressure * unit_normal(0) * a_plus;
      break;
    }

    flux_dot_n(ConservativeVariables::TOTAL_ENERGY_DENSITY) = 
      flux_dot_n(ConservativeVariables::MASS_DENSITY) / conservative_state[ConservativeVariables::MASS_DENSITY] *
      ( conservative_state[ConservativeVariables::TOTAL_ENERGY_DENSITY] + pressure ) - 1. / 2. * pressure * b_plus;

    flux_dot_n *= d_area; // rescale final result

    // Compute maximum characteristic speed

    // sound speed
    const double sound = euler::speedOfSound(species_, conservative_state[ConservativeVariables::MASS_DENSITY], pressure);
    // fluid speed
    const double speed = std::fabs(normal_velocity);
    // max characteristic speed = fluid speed + sound speed
    return speed + sound;
  }

} // namespace mfpic
