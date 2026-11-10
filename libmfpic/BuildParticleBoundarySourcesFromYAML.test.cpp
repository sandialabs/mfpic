#include <libmfpic/BuildParticleBoundarySourcesFromYAML.hpp>
#include <libmfpic/Constants.hpp>
#include <libmfpic/MaxwellianParticleBoundarySource.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace {

using namespace mfpic;
using namespace testing;

const Species electron_species{.charge = -constants::elementary_charge, .mass = constants::electron_mass, .name = "electron"};
const Species proton_species{.charge = constants::elementary_charge, .mass = constants::proton_mass, .name = "proton"};

const std::unordered_map<std::string, Species> species_map{
  {"electron", electron_species},
  {"proton", proton_species}};

constexpr int velocity_dims = 3;

std::shared_ptr<mfem::Mesh> makeUnitMesh(int mesh_dimension) {
  switch (mesh_dimension) {
    case 1: return std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(2));
    case 2: return std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian2D(2, 2, mfem::Element::QUADRILATERAL));
    default: return std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::HEXAHEDRON));
  }
}

std::vector<std::unique_ptr<ParticleBoundarySource>> buildFromString(const std::string& yaml_string, int mesh_dimension = 1) {
  const YAML::Node node = YAML::Load(yaml_string);
  return buildParticleBoundarySourcesFromYAML(node["Boundary Sources"], species_map, makeUnitMesh(mesh_dimension), velocity_dims);
}

const MaxwellianParticleBoundarySource& asMaxwellian(const std::unique_ptr<ParticleBoundarySource>& boundary_source) {
  const auto* maxwellian = dynamic_cast<const MaxwellianParticleBoundarySource*>(boundary_source.get());
  if (maxwellian == nullptr) {
    throw std::runtime_error("Boundary source is not Maxwellian");
  }
  return *maxwellian;
}

void expectThrowsWithMessage(const std::string& yaml_string, const std::string& message) {
  EXPECT_THAT(
    [&]() {buildFromString(yaml_string);},
    ThrowsMessage<std::logic_error>(HasSubstr(message))
  );
}

TEST(ParticleBoundarySourcesFactory, MissingBoundarySourcesGivesEmptyList) {
  EXPECT_TRUE(buildFromString("Other: 1\n").empty());
}

TEST(ParticleBoundarySourcesFactory, NullBoundarySourcesGivesEmptyList) {
  EXPECT_TRUE(buildFromString("Boundary Sources:\n").empty());
}

TEST(ParticleBoundarySourcesFactory, OneSourceWithTwoSpeciesGivesTwoMaxwellianParameters) {
  const auto boundary_sources = buildFromString(
    "Boundary Sources:\n"
    "  - Species: [electron, proton]\n"
    "    Side: right\n"
    "    Number of Macroparticles per Species per Timestep: 7\n"
    "    Maxwellian:\n"
    "      Flux Density: 2.5e20\n"
    "      Temperature: 300.\n");

  ASSERT_EQ(boundary_sources.size(), 2u);
  const std::vector<Species> expected_species{electron_species, proton_species};
  for (int i = 0; i < 2; i++) {
    const auto* maxwellian = dynamic_cast<const MaxwellianParticleBoundarySource*>(boundary_sources[i].get());
    ASSERT_NE(maxwellian, nullptr);
    EXPECT_EQ(maxwellian->getSpecies(), expected_species[i]);
    EXPECT_EQ(maxwellian->getBoundaryAttribute(), 2);
    EXPECT_EQ(maxwellian->getNumParticlesPerTimestep(), 7);
    EXPECT_DOUBLE_EQ(maxwellian->getFluxDensity(), 2.5e20);
    EXPECT_DOUBLE_EQ(maxwellian->getTemperature(), 300.);
  }
}

TEST(ParticleBoundarySourcesFactory, IntegerSideIsBoundaryAttribute) {
  const auto boundary_sources = buildFromString(
    "Boundary Sources:\n"
    "  - Species: [electron]\n"
    "    Side: 3\n"
    "    Number of Macroparticles per Species per Timestep: 1\n"
    "    Maxwellian:\n"
    "      Flux Density: 1.\n"
    "      Temperature: 1.\n",
    2);

  EXPECT_EQ(asMaxwellian(boundary_sources[0]).getBoundaryAttribute(), 3);
}

TEST(ParticleBoundarySourcesFactory, InvalidSideThrows) {
  expectThrowsWithMessage(
    "Boundary Sources:\n"
    "  - Species: [electron]\n"
    "    Side: nowhere\n"
    "    Number of Macroparticles per Species per Timestep: 1\n"
    "    Maxwellian:\n"
    "      Flux Density: 1.\n"
    "      Temperature: 1.\n",
    "Side must be an inline mesh side name or a positive integer boundary attribute!");
}

TEST(ParticleBoundarySourcesFactory, MissingSideThrows) {
  expectThrowsWithMessage(
    "Boundary Sources:\n"
    "  - Species: [electron]\n"
    "    Number of Macroparticles per Species per Timestep: 1\n"
    "    Maxwellian:\n"
    "      Flux Density: 1.\n"
    "      Temperature: 1.\n",
    "Side is required!");
}

TEST(ParticleBoundarySourcesFactory, SpeciesMustBeASequence) {
  expectThrowsWithMessage(
    "Boundary Sources:\n"
    "  - Species: electron\n"
    "    Side: left\n"
    "    Number of Macroparticles per Species per Timestep: 1\n"
    "    Maxwellian:\n"
    "      Flux Density: 1.\n"
    "      Temperature: 1.\n",
    "Species is not a sequence!");
}

TEST(ParticleBoundarySourcesFactory, SpeciesMustExistInSpeciesMap) {
  expectThrowsWithMessage(
    "Boundary Sources:\n"
    "  - Species: [nonexistent]\n"
    "    Side: left\n"
    "    Number of Macroparticles per Species per Timestep: 1\n"
    "    Maxwellian:\n"
    "      Flux Density: 1.\n"
    "      Temperature: 1.\n",
    "Species \"nonexistent\" was not created in the Species block!");
}

TEST(ParticleBoundarySourcesFactory, NumberOfMacroparticlesIsRequired) {
  expectThrowsWithMessage(
    "Boundary Sources:\n"
    "  - Species: [electron]\n"
    "    Side: left\n"
    "    Maxwellian:\n"
    "      Flux Density: 1.\n"
    "      Temperature: 1.\n",
    "Number of Macroparticles per Species per Timestep is required for a boundary source!");
}

TEST(ParticleBoundarySourcesFactory, NumberOfMacroparticlesMustBePositive) {
  expectThrowsWithMessage(
    "Boundary Sources:\n"
    "  - Species: [electron]\n"
    "    Side: left\n"
    "    Number of Macroparticles per Species per Timestep: 0\n"
    "    Maxwellian:\n"
    "      Flux Density: 1.\n"
    "      Temperature: 1.\n",
    "Number of Macroparticles per Species per Timestep is nonpositive!");
}

TEST(ParticleBoundarySourcesFactory, FluxDensityMustBePositive) {
  expectThrowsWithMessage(
    "Boundary Sources:\n"
    "  - Species: [electron]\n"
    "    Side: left\n"
    "    Number of Macroparticles per Species per Timestep: 1\n"
    "    Maxwellian:\n"
    "      Flux Density: 0.\n"
    "      Temperature: 1.\n",
    "Flux Density is nonpositive!");
}

TEST(ParticleBoundarySourcesFactory, TemperatureIsRequired) {
  expectThrowsWithMessage(
    "Boundary Sources:\n"
    "  - Species: [electron]\n"
    "    Side: left\n"
    "    Number of Macroparticles per Species per Timestep: 1\n"
    "    Maxwellian:\n"
    "      Flux Density: 1.\n",
    "Temperature is required!");
}

TEST(ParticleBoundarySourcesFactory, TemperatureMustBePositive) {
  expectThrowsWithMessage(
    "Boundary Sources:\n"
    "  - Species: [electron]\n"
    "    Side: left\n"
    "    Number of Macroparticles per Species per Timestep: 1\n"
    "    Maxwellian:\n"
    "      Flux Density: 1.\n"
    "      Temperature: -1.\n",
    "Temperature is nonpositive!");
}

TEST(ParticleBoundarySourcesFactory, BoundarySourceTypeIsRequired) {
  expectThrowsWithMessage(
    "Boundary Sources:\n"
    "  - Species: [electron]\n"
    "    Side: left\n"
    "    Number of Macroparticles per Species per Timestep: 1\n",
    "Missing a boundary source type.");
}

} // namespace
