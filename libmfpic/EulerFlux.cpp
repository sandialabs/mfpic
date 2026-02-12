#include <libmfpic/Euler.hpp>
#include <libmfpic/EulerFlux.hpp>

namespace mfpic {

  double EulerFlux::ComputeFlux(const mfem::Vector &conservative_state,
                                mfem::ElementTransformation &,
                                mfem::DenseMatrix &flux) const
  {
    using namespace euler;
    const double pressure = getPressureFromConservativeState(conservative_state, species_);
    const mfem::Vector momentum = getMomentumDensityFromConservativeState(conservative_state);
    const mfem::Vector velocity = getBulkVelocityFromConservativeState(conservative_state);

    // Check whether the solution is physical only in debug mode
    MFEM_ASSERT(conservative_state[ConservativeVariables::MASS_DENSITY] >= 0, "Negative Density");
    MFEM_ASSERT(pressure >= 0, "Negative Pressure");
    MFEM_ASSERT(conservative_state[ConservativeVariables::TOTAL_ENERGY_DENSITY] >= 0, "Negative Total Energy Density");

     for (int d = 0; d < dim; d++)
     {
        flux(ConservativeVariables::MASS_DENSITY, d) = momentum(d);

        flux(ConservativeVariables::X_MOMENTUM_DENSITY, d) = 
          momentum(0) * velocity(d) + double(d==0) * pressure;
        flux(ConservativeVariables::Y_MOMENTUM_DENSITY, d) = 
          momentum(1) * velocity(d) + double(d==1) * pressure;
        flux(ConservativeVariables::Z_MOMENTUM_DENSITY, d) = 
          momentum(2) * velocity(d) + double(d==2) * pressure;

        flux(ConservativeVariables::TOTAL_ENERGY_DENSITY, d) = 
          velocity(d) * (pressure + conservative_state[ConservativeVariables::TOTAL_ENERGY_DENSITY]);
     }

    // Compute maximum characteristic speed

    // sound speed
    const double sound = euler::speedOfSound(species_, conservative_state[ConservativeVariables::MASS_DENSITY], pressure);
    // fluid speed
    const double speed = velocity.Norml2();
    // max characteristic speed = fluid speed + sound speed
    return speed + sound;
  }


  double EulerFlux::ComputeFluxDotN(const mfem::Vector &conservative_state, 
                                    const mfem::Vector &normal,
                                    mfem::FaceElementTransformations &,
                                    mfem::Vector &flux_dot_n) const
  {
    using namespace euler;
    const double pressure = getPressureFromConservativeState(conservative_state, species_);
    const mfem::Vector momentum = getMomentumDensityFromConservativeState(conservative_state);
    const mfem::Vector velocity = getBulkVelocityFromConservativeState(conservative_state);

    // Check whether the solution is physical only in debug mode
    MFEM_ASSERT(conservative_state[ConservativeVariables::MASS_DENSITY] >= 0, "Negative Density");
    MFEM_ASSERT(pressure >= 0, "Negative Pressure");
    MFEM_ASSERT(conservative_state[ConservativeVariables::TOTAL_ENERGY_DENSITY] >= 0, "Negative Total Energy Density");

    double normal_velocity(0.);
    for (int d = 0; d < dim; d++)
      normal_velocity += velocity(d) * normal(d);

    flux_dot_n(ConservativeVariables::MASS_DENSITY) = conservative_state[ConservativeVariables::MASS_DENSITY] * normal_velocity;

    flux_dot_n(ConservativeVariables::X_MOMENTUM_DENSITY) = normal_velocity * momentum(0);
    flux_dot_n(ConservativeVariables::Y_MOMENTUM_DENSITY) = normal_velocity * momentum(1);
    flux_dot_n(ConservativeVariables::Z_MOMENTUM_DENSITY) = normal_velocity * momentum(2);

    switch (dim) {
    case 3:
      flux_dot_n(ConservativeVariables::Z_MOMENTUM_DENSITY) += pressure * normal(2);
      [[fallthrough]];
    case 2:
      flux_dot_n(ConservativeVariables::Y_MOMENTUM_DENSITY) += pressure * normal(1);
      [[fallthrough]];
    case 1:
      flux_dot_n(ConservativeVariables::X_MOMENTUM_DENSITY) += pressure * normal(0);
      break;
    }

    flux_dot_n(ConservativeVariables::TOTAL_ENERGY_DENSITY) = 
      normal_velocity * (conservative_state[ConservativeVariables::TOTAL_ENERGY_DENSITY] + pressure);

    // Compute maximum characteristic speed

    // sound speed
    const double sound = euler::speedOfSound(species_, conservative_state[ConservativeVariables::MASS_DENSITY], pressure);
    // fluid speed
    const double speed = std::fabs(normal_velocity) / std::sqrt(normal*normal);
    // max characteristic speed = fluid speed + sound speed
    return speed + sound;
  }

} // namespace mfpic
