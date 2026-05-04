#pragma once

#include <libmfpic/Species.hpp>

#include <mfem/mfem.hpp>
#include <yaml-cpp/yaml.h>

namespace mfpic {

/**
 * @brief parameters that define the state of a species at a single location.
 */
struct SourceStateParameters {
  double number_density;
  mfem::Vector bulk_velocity{0., 0., 0.};
  double temperature = 0;
  double kappa = -1;
};

/**
 * @brief parameters that define a species over the whole domain, can either be used to load a species as a source or as initial
 *  conditions.
 */
struct SourceParameters {
  Species species;
  int num_particles = 0;

  SourceParameters() {};
  SourceParameters(
    const Species& species_in,
    const int num_particles_in = 0)
    : species(species_in)
    , num_particles(num_particles_in){};

  virtual ~SourceParameters() = default;

  /**
   * @brief Evaluate the basic state parameters at a point in space.
   *
   * @param[in] x Point in space at which to evaluate.
   *
   * @returns Basic state parameters at @a x.
   */
  virtual SourceStateParameters sourceStateParametersAtPoint(const mfem::Vector& x) const = 0;

  /**
   * @brief Get an mfem::VectorFunctionCoefficient that represents an Euler fluid with the parameters in this object
   * 
   * @return mfem::VectorFunctionCoefficient
   */
  mfem::VectorFunctionCoefficient getEulerVectorCoefficient() const;
};

/**
 * @brief parameters defining a species that is constant in space
 */
struct ConstantSourceParameters : public SourceParameters {
  SourceStateParameters constant_state;

  ConstantSourceParameters(
    const Species& species,
    const SourceStateParameters& state_parameters,
    const int num_particles = 0);

  ConstantSourceParameters(
    const Species& species,
    const double number_density,
    const double temperature,
    const mfem::Vector& bulk_velocity=mfem::Vector({0., 0., 0.}),
    const double kappa=-1);

  /**
   * @brief Give the basic state parameters that represent a uniform distribution in space.
   *
   * @returns The constant state parameters.
   */
  virtual SourceStateParameters sourceStateParametersAtPoint(const mfem::Vector&) const override;
};

/**
 * @brief parameters defining a species that is piecewise constant with one discontinuity 
 *  if x < discontinuity_location then species is given by left state otherwise species is right state
 */
struct SodSourceParameters : public SourceParameters {
  double discontinuity_location;
  SourceStateParameters left_state;
  SourceStateParameters right_state;

  SodSourceParameters(
    const Species& species,
    const double discontinuity_location,
    const SourceStateParameters& left_state_parameters,
    const SourceStateParameters& right_state_parameters,
    const int num_particles = 0);

  /**
   * @brief Evaluate the basic state parameters that represent a piecewise uniform distribution.
   *
   * @param[in] x Point in space at which to evaluate.
   *
   * @returns Basic state parameters at @a x.
   */
  virtual SourceStateParameters sourceStateParametersAtPoint(const mfem::Vector& x) const override;
};

/**
 * @brief parameters defining a species that is a gaussian in space, i.e. the number density is given by
 *  n = h_n exp(-0.5 (x - mu)^2 / sigma) + o_n
 *  where h_n is the height of the gaussian, o_n is the offset, mu is the center or mean of the gaussian, and sigma is the 
 *  standard deviation.
 *  The formulas for the velocity and pressure are the same except that they use their own heights and offsets
 */
struct GaussianSourceParameters: public SourceParameters {
  mfem::Vector center;
  double standard_deviation;

  // temperature is ignored in these states
  SourceStateParameters offsets;
  SourceStateParameters heights;

  double pressure_offset;
  double pressure_height;

  GaussianSourceParameters(
    const Species& species,
    const mfem::Vector& center,
    const double standard_deviation,
    const SourceStateParameters& offsets,
    const SourceStateParameters& heights,
    const double pressure_offset,
    const double pressure_height,
    const int num_particles = 0);

  /**
   * @brief Evaluate the basic state parameters that represent a distribution that is Gaussian in space.
   *
   * @param[in] x Point in space at which to evaluate.
   *
   * @returns Basic state parameters at @a x.
   */
  virtual SourceStateParameters sourceStateParametersAtPoint(const mfem::Vector& x) const override;
};

/**
 * @brief parameters defining a species that is nearly constant in space but 
 * perturbed in a periodic sense. The wavenumber \f$k\f$ and perturbation amplitude \f$\eps\f$ are supplied
 * to give a form \f$g = C (1 + \eps cos(k x)) \f$.
 */
struct PeriodicPerturbationSourceParameters : public SourceParameters {
  mfem::Vector wavevector;
  SourceStateParameters base_values;
  SourceStateParameters perturbations;

  PeriodicPerturbationSourceParameters(
    const Species& species,
    const mfem::Vector& wavevector,
    const SourceStateParameters& base_values,
    const SourceStateParameters& perturbations,
    const int num_particles = 0);

  /**
   * @brief Evaluate the basic state parameters that represent a distribution that is periodically perturbed in space.
   *
   * @param[in] x Point in space at which to evaluate.
   *
   * @returns Basic state parameters at @a x.
   */
  virtual SourceStateParameters sourceStateParametersAtPoint(const mfem::Vector& x) const override;
};

/**
 * @brief build SourceStateParameters that defines a state from YAML
 * 
 * @param state_node - the YAML defining the state
 * @return SourceStateParameters - the parameters defining the state for a source or initial condition
 */
SourceStateParameters buildSourceStateParametersFromYAML(const YAML::Node& state_node);

/**
 * @brief build a list of parameters that define sources or initial conditions from YAML
 * 
 * @param sources_node - the YAML defining the sources or initial conditions
 * @param species_map - a map of species name to species definition
 * @return std::vector<std::unique_ptr<SourceParameters>> - a list of parameters defining the sources or initial conditions
 */
std::vector<std::unique_ptr<SourceParameters>> buildListOfSourceParametersFromYAML(
  const YAML::Node& sources_node,
  const std::unordered_map<std::string, Species>& species_map);

/**
 * @brief build a list of coefficients from a vector of \ref SourceParameters
 *
 * @param list_of_parameters a list of parameters defining the sources or initial conditions
 * @returns std::vector<std::pair<Species, std::unique_ptr<mfem::VectorCoefficient>>> - a list of \ref Species and vector coefficient pairs 
 *
 * @note Appropriate for Euler fluids. The \ref SourceParameters passed in must support \ref getEulerVectorCoefficient .
 */
std::vector<std::pair<Species, std::unique_ptr<mfem::VectorCoefficient>>> buildListOfSpeciesAndEulerSourceCoefficients(
  const std::vector<std::unique_ptr<SourceParameters>>& list_of_parameters);

}
