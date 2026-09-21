#include <libmfpic/BuildParticlesFromYaml.hpp>
#include <libmfpic/RandomNumberGenerator.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <mfem/mfem.hpp>
#include <yaml-cpp/yaml.h>

#include <random>

namespace {

using namespace mfpic;
using namespace testing;

constexpr int num_elements = 10;
const std::shared_ptr<mfem::Mesh> simple_mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elements));
const std::unordered_map<std::string, Species> species_map_with_electron = {
  {"electron", Species{.charge = 1.6e-19, .mass = 9.11e-31}},
};

TEST(BuildParticlesFromYaml, NumberOfParticlesLoadedIsAsExpected) {
  std::unordered_map<std::string, Species> species_map = species_map_with_electron;
  species_map["fake_proton"] = Species{};
  const std::string yaml = R"(
Initial Conditions:
  - Species: [electron, fake_proton]
    Number of Macroparticles per Species: 2
    Constant:
      Temperature: 1.0
      Number Density: 1.0
  - Species: [electron]
    Number of Macroparticles per Species: 4
    Constant:
      Temperature: 1.0
      Number Density: 1.0
  )";
  const YAML::Node node = YAML::Load(yaml);

  RandomNumberGenerator generator;
  ParticleContainer particles = buildParticlesFromYaml(
    node["Initial Conditions"],
    species_map,
    generator,
    simple_mesh,
    3
  );
  EXPECT_EQ(particles.numParticles(), 8);
}

TEST(BuildParticlesFromYaml, NumberOfNonZeroVelocityDimensionsIsCorrect) {
  std::unordered_map<std::string, Species> species_map = species_map_with_electron;
  species_map["fake_proton"] = Species{};
  const std::string yaml = R"(
Initial Conditions:
  - Species: [electron, fake_proton]
    Number of Macroparticles per Species: 2
    Constant:
      Temperature: 1.0
      Number Density: 1.0
  - Species: [electron]
    Number of Macroparticles per Species: 4
    Constant:
      Temperature: 1.0
      Number Density: 1.0
  )";
  const YAML::Node node = YAML::Load(yaml);

  RandomNumberGenerator generator;
  ParticleContainer particles = buildParticlesFromYaml(
    node["Initial Conditions"],
    species_map,
    generator,
    simple_mesh,
    1
  );
  EXPECT_EQ(particles.numParticles(), 8);
  for (Particle particle : particles)
  {
    EXPECT_EQ(particle.velocity(1), 0.0);
    EXPECT_EQ(particle.velocity(2), 0.0);
  }

  particles = buildParticlesFromYaml(
    node["Initial Conditions"],
    species_map,
    generator,
    simple_mesh,
    2
  );

  EXPECT_EQ(particles.numParticles(), 8);
  for (Particle particle : particles)
  {
    EXPECT_EQ(particle.velocity(2), 0.0);
  }
}

} // namespace
