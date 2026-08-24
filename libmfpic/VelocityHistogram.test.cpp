#include <libmfpic/Constants.hpp>
#include <libmfpic/ParticleContainer.hpp>
#include <libmfpic/ParticleOperations.hpp>

#include <gtest/gtest.h>
#include <mfem/mesh/element.hpp>
#include <unordered_map>

namespace {

using namespace mfpic;

TEST(ParticleOperations, VelocityHistogramIntegratesToNumberDensity)
{
  Species species{.charge = -constants::elementary_charge, .mass = constants::electron_mass};
  std::mt19937 generator;

  const int num_elems = 5;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems,1.0);

  ParticleContainer particles;

  Species electron {.mass = constants::electron_mass};

  particles.addParticle({
      .velocity = mfem::Vector({-1.0, 0.0, 0.0}),
      .species = electron,
      .weight = 1.0,
      .is_alive = true
  });

  particles.addParticle({
      .velocity = mfem::Vector({0.0, 0.0, 0.0}),
      .species = electron,
      .weight = 2.0,
      .is_alive = true
  });

  particles.addParticle({
      .velocity = mfem::Vector({1.0, 0.0, 0.0}),
      .species = electron,
      .weight = 1.0,
      .is_alive = true
  });

  particles.addParticle({
      .velocity = mfem::Vector({100.0, 0.0, 0.0}),
      .species = electron,
      .weight = 100.0,
      .is_alive = false
  });

  constexpr int nbins = 4;

  VelocityHistogram hist;
  hist.buildVelocityHistogram(particles, mesh, nbins);

  EXPECT_EQ(hist.getNumBins(), nbins);
  EXPECT_EQ(hist.getNumEdges(), nbins + 1);

  EXPECT_DOUBLE_EQ(hist.getVMin(), -1.0);
  EXPECT_DOUBLE_EQ(hist.getVMax(), 1.0);

  EXPECT_DOUBLE_EQ(hist.getDV(), 0.5);

  double integral = 0.0;
  for (double value : hist.getBinValues())
      integral += value * hist.getDV();

  EXPECT_NEAR(integral, 4.0, 1e-12);
}

TEST(ParticleOperations, UpdatesParticleDistributionFunctionValueFromHistogram)
{
  ParticleContainer particles;

  Species electron {.mass = constants::electron_mass};

  particles.addParticle({
      .velocity = mfem::Vector({-0.75, 0.0, 0.0}),
      .species = electron,
      .weight = 1.0,
      .is_alive = true
  });

  particles.addParticle({
      .velocity = mfem::Vector({0.75, 0.0, 0.0}),
      .species = electron,
      .weight = 1.0,
      .is_alive = true
  });

  const int num_elems = 5;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems,1.0);

  VelocityHistogram hist;
  hist.buildVelocityHistogram(particles,mesh,2);

  ParticleOperations::updateParticleDistributionFunctionValue(
      particles,
      hist);

  int index = 0;
  for (const Particle& particle : particles)
  {
      if (!particle.is_alive)
          continue;

      EXPECT_DOUBLE_EQ(
          particle.particle_distribution_function_value,
          hist.evaluate(particle.velocity(0)));

      index++;
  }

  EXPECT_EQ(index, 2);
}

} // namespace
