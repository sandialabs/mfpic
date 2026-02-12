#include <libmfpic/BuildSpeciesMapFromYaml.hpp>

#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

namespace {

using namespace mfpic;

TEST(BuildSpeciesMapFromYaml, ThrowWhenSpeciesIsNotAMap) {
  const std::string not_a_map = "[3.0, 3.0, 4.0]";
  YAML::Node node = YAML::Load(not_a_map);
  EXPECT_ANY_THROW(buildSpeciesMapFromYaml(node));
}

TEST(BuildSpeciesMapFromYaml, ChargeOverMassComputedCorrectlyWhenNotExplicitlyProvided) {
  const std::string yaml = "electron: {Charge: 1.0, Mass: 2.0}";
  YAML::Node node = YAML::Load(yaml);

  std::unordered_map<std::string, Species> species_map = buildSpeciesMapFromYaml(node);

  EXPECT_EQ(species_map.size(), 1);
  for (auto [name, species] : species_map) {
    EXPECT_EQ(name, "electron");
    EXPECT_DOUBLE_EQ(species.charge, 1.0);
    EXPECT_DOUBLE_EQ(species.mass, 2.0);
    EXPECT_DOUBLE_EQ(species.charge_over_mass, 0.5);
  }
}

TEST(BuildSpeciesMapFromYaml, ChargeOverMassGotWhenExplicitlyProvided) {
  const std::string yaml = "electron: {Charge: 1.0, Mass: 2.0, Charge Over Mass: 3.0}";
  YAML::Node node = YAML::Load(yaml);

  std::unordered_map<std::string, Species> species_map = buildSpeciesMapFromYaml(node);

  EXPECT_EQ(species_map.size(), 1);
  for (auto [name, species] : species_map) {
    EXPECT_EQ(name, "electron");
    EXPECT_DOUBLE_EQ(species.charge, 1.0);
    EXPECT_DOUBLE_EQ(species.mass, 2.0);
    EXPECT_DOUBLE_EQ(species.charge_over_mass, 3.0);
  }
}

TEST(BuildSpeciesMapFromYaml, MultipleSpeciesAreGrabbed) {
  const std::string yaml = R"(
electron:
  Charge: 1.0
  Mass: 2.0
proton:
  Charge: 3.0
  Mass: 4.0
  )";
  YAML::Node node = YAML::Load(yaml);

  std::unordered_map<std::string, Species> species_map = buildSpeciesMapFromYaml(node);

  EXPECT_EQ(species_map.size(), 2);
  const Species electron = species_map.at("electron");
  EXPECT_DOUBLE_EQ(electron.charge, 1.0);
  EXPECT_DOUBLE_EQ(electron.mass, 2.0);
  const Species proton = species_map.at("proton");
  EXPECT_DOUBLE_EQ(proton.charge, 3.0);
  EXPECT_DOUBLE_EQ(proton.mass, 4.0);
}

TEST(BuildSpeciesMapFromYaml, SpecificHeatRatioDefaultsToFiveThirds) {
  const std::string yaml = R"(
electron:
  Charge: 1.0
  Mass: 2.0
)";

  YAML::Node node = YAML::Load(yaml);
  std::unordered_map<std::string, Species> species_map = buildSpeciesMapFromYaml(node);

  EXPECT_EQ(species_map.size(), 1);
  for (auto [name, species] : species_map) {
    EXPECT_DOUBLE_EQ(species.specific_heat_ratio, 5. / 3.);
  }
}

TEST(BuildSpeciesMapFromYaml, SpecificHeatRatiosIsPickedUpWhenSpecified) {
  constexpr double specific_heat_ratio = 1.4;
  const std::string yaml = R"(
air:
  Charge: 1.0
  Mass: 2.0
  Specific Heat Ratio: )" + std::to_string(specific_heat_ratio) + R"(
)";

  YAML::Node node = YAML::Load(yaml);
  std::unordered_map<std::string, Species> species_map = buildSpeciesMapFromYaml(node);

  EXPECT_EQ(species_map.size(), 1);
  for (auto [name, species] : species_map) {
    EXPECT_DOUBLE_EQ(species.specific_heat_ratio, specific_heat_ratio);
  }
}

} // namespace
