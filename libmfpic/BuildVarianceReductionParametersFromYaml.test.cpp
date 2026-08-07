#include <libmfpic/BuildVarianceReductionParametersFromYaml.hpp>

#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

namespace {

using namespace mfpic;

TEST(BuildVarianceReductionParametersFromYaml, VarianceReductionParametersSetCorrectly) {
  const std::string yaml = R"(
Strategy: Euler Fluid
Specified LF Moments:
  Reference Number Density: 9320.59320
  Reference Bulk Velocity: [53920.183,620.1988,28.190]
  Reference Temperature: 53290.201
Use VR E Field: true
Apply Cell Limiting: false
  )";
  const YAML::Node node = YAML::Load(yaml);

  auto params= buildVarianceReductionParametersFromYAML(node);

  EXPECT_EQ(params.strategy,VarianceReductionParameters::Strategy::EulerFluid);
  EXPECT_EQ(params.limit_variance_reduction,false);
  EXPECT_EQ(params.use_variance_reduced_electric_field,true);
  EXPECT_EQ(params.reference_number_density,9320.59320);
  EXPECT_EQ(params.reference_bulk_velocity[0],53920.183);
  EXPECT_EQ(params.reference_bulk_velocity[1],620.1988);
  EXPECT_EQ(params.reference_bulk_velocity[2],28.190);
  EXPECT_EQ(params.reference_temperature,53290.201);
}

}
