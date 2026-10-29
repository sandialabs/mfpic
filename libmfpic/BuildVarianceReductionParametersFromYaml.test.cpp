#include <libmfpic/BuildVarianceReductionParametersFromYaml.hpp>

#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

namespace {

using namespace mfpic;

TEST(BuildVarianceReductionParametersFromYaml, VarianceReductionParametersSetCorrectly) {
  const std::string yaml = R"(
Strategy: Uniform Maxwellian
Uniform Maxwellian Moments:
  Number Density: 9320.59320
  Bulk Velocity: [53920.183,620.1988,28.190]
  Temperature: 53290.201
Use VR E Field: true
Apply Cell Limiting: true
  )";
  const YAML::Node node = YAML::Load(yaml);

  auto params= buildVarianceReductionParametersFromYAML(node);

  EXPECT_EQ(params.strategy,VarianceReductionParameters::Strategy::UniformMaxwellian);
  EXPECT_EQ(params.limit_variance_reduction,true);
  EXPECT_EQ(params.use_variance_reduced_electric_field,true);
  EXPECT_EQ(params.reference_number_density,9320.59320);
  EXPECT_EQ(params.reference_bulk_velocity[0],53920.183);
  EXPECT_EQ(params.reference_bulk_velocity[1],620.1988);
  EXPECT_EQ(params.reference_bulk_velocity[2],28.190);
  EXPECT_EQ(params.reference_temperature,53290.201);
  EXPECT_TRUE(params.reference_species.empty());
}

TEST(BuildVarianceReductionParametersFromYaml, UniformMaxwellianSpeciesAreParsed) {
  const std::string yaml = R"(
Strategy: Uniform Maxwellian
Uniform Maxwellian Moments:
  Number Density: 1.0
  Bulk Velocity: [0.0, 0.0, 0.0]
  Temperature: 1.0
  Species: [electron, proton]
  )";
  const YAML::Node node = YAML::Load(yaml);

  auto params = buildVarianceReductionParametersFromYAML(node);

  EXPECT_EQ(params.reference_species, (std::vector<std::string>{"electron", "proton"}));
}

TEST(BuildVarianceReductionParametersFromYaml, DefaultsAreEulerFluidWithoutCellLimiting) {
  const YAML::Node node = YAML::Load("Use VR E Field: true");

  auto params = buildVarianceReductionParametersFromYAML(node);

  EXPECT_EQ(params.strategy, VarianceReductionParameters::Strategy::EulerFluid);
  EXPECT_EQ(params.limit_variance_reduction, false);
}

TEST(BuildVarianceReductionParametersFromYaml, NoneStrategyIsParsed) {
  const YAML::Node node = YAML::Load("Strategy: none");

  EXPECT_EQ(buildVarianceReductionParametersFromYAML(node).strategy, VarianceReductionParameters::Strategy::None);
}

TEST(BuildVarianceReductionParametersFromYaml, UnknownStrategyThrows) {
  const YAML::Node node = YAML::Load("Strategy: Spatially Averaged");

  EXPECT_THROW(buildVarianceReductionParametersFromYAML(node), std::logic_error);
}

TEST(BuildVarianceReductionParametersFromYaml, UniformMaxwellianWithoutMomentsThrows) {
  const YAML::Node node = YAML::Load("Strategy: Uniform Maxwellian");

  EXPECT_THROW(buildVarianceReductionParametersFromYAML(node), std::logic_error);
}

TEST(BuildVarianceReductionParametersFromYaml, UniformMaxwellianMomentsWithEulerFluidThrows) {
  const std::string yaml = R"(
Strategy: Euler Fluid
Uniform Maxwellian Moments:
  Number Density: 1.0
  Bulk Velocity: [0.0, 0.0, 0.0]
  Temperature: 1.0
  )";
  const YAML::Node node = YAML::Load(yaml);

  EXPECT_THROW(buildVarianceReductionParametersFromYAML(node), std::logic_error);
}

TEST(BuildVarianceReductionParametersFromYaml, UniformMaxwellianBulkVelocityMustHaveLengthThree) {
  const std::string yaml = R"(
Strategy: Uniform Maxwellian
Uniform Maxwellian Moments:
  Number Density: 1.0
  Bulk Velocity: [0.0]
  Temperature: 1.0
  )";
  const YAML::Node node = YAML::Load(yaml);

  EXPECT_THROW(buildVarianceReductionParametersFromYAML(node), std::logic_error);
}

TEST(BuildVarianceReductionParametersFromYaml, SpecifiedLFMomentsIsRejected) {
  const std::string yaml = R"(
Strategy: Euler Fluid
Specified LF Moments:
  Reference Number Density: 1.0
  Reference Bulk Velocity: [0.0, 0.0, 0.0]
  Reference Temperature: 1.0
  )";
  const YAML::Node node = YAML::Load(yaml);

  EXPECT_THROW(buildVarianceReductionParametersFromYAML(node), std::logic_error);
}

}
