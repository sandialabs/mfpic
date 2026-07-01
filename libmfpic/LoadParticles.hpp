#pragma once

#include <libmfpic/Euler.hpp>
#include <libmfpic/GenerateKappaVelocity.hpp>
#include <libmfpic/GenerateMaxwellianVelocity.hpp>
#include <libmfpic/MeshDistribution.hpp>
#include <libmfpic/MeshUtilities.hpp>
#include <libmfpic/ParticleContainer.hpp>
#include <libmfpic/RandomNumberGenerator.hpp>
#include <libmfpic/SourcesFactory.hpp>
#include <libmfpic/Species.hpp>

#include <mfem/mfem.hpp>

namespace mfpic {

/**
 * @brief Load a requested number of particles according to a parametrized distribution.
 *
 * @param[in]     source_parameters          Parameters for the particle distribution.
 * @param[in,out] generator                  Random number generator.
 * @param[in]     mesh                       Mesh in which to create particles.
 *
 * @returns Container of created particles.
 */
ParticleContainer loadParticles(
  const SourceParameters& source_parameters,
  RandomNumberGenerator& generator,
  std::shared_ptr<mfem::Mesh> mesh
);

} // namespace mfpic
