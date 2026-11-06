#include <libmfpic/Constants.hpp>
#include <libmfpic/ParseHelpers.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace {

using namespace mfpic;
using namespace testing;

const Species proton_species{.charge = constants::elementary_charge, .mass = constants::proton_mass, .name = "proton"};
const std::unordered_map<std::string, Species> species_map{{"proton", proton_species}};

template <typename Function>
void expectThrowsWithMessage(Function function, const std::string& message) {
  EXPECT_THAT(function, ThrowsMessage<std::logic_error>(HasSubstr(message)));
}

TEST(ParseHelpers, BoundaryAttributeFromInlineSideNameDependsOnMeshDimension) {
  const YAML::Node node = YAML::Load("Side: left");
  EXPECT_EQ(parseBoundaryAttribute(node, 1), 1);
  EXPECT_EQ(parseBoundaryAttribute(node, 2), 4);
  EXPECT_EQ(parseBoundaryAttribute(node, 3), 5);
}

TEST(ParseHelpers, BoundaryAttributeFromInteger) {
  EXPECT_EQ(parseBoundaryAttribute(YAML::Load("Side: 7"), 1), 7);
}

TEST(ParseHelpers, MissingSideThrows) {
  expectThrowsWithMessage([]() { parseBoundaryAttribute(YAML::Load("Other: 1"), 1); }, "Side is required!");
}

TEST(ParseHelpers, InvalidSideThrows) {
  const std::string message = "Side must be an inline mesh side name or a positive integer boundary attribute!";
  expectThrowsWithMessage([]() { parseBoundaryAttribute(YAML::Load("Side: nowhere"), 1); }, message);
  expectThrowsWithMessage([]() { parseBoundaryAttribute(YAML::Load("Side: 0"), 1); }, message);
  expectThrowsWithMessage([]() { parseBoundaryAttribute(YAML::Load("Side: -2"), 1); }, message);
}

TEST(ParseHelpers, SpeciesNameMustExist) {
  EXPECT_EQ(parseSpeciesName(YAML::Load("proton"), species_map), "proton");
  expectThrowsWithMessage(
    []() { parseSpeciesName(YAML::Load("electron"), species_map); },
    "Species was not created in the Species block!");
}

TEST(ParseHelpers, SpeciesNamesMustBeASequenceOfExistingSpecies) {
  EXPECT_THAT(parseSpeciesNames(YAML::Load("[proton, proton]"), species_map), ElementsAre("proton", "proton"));
  expectThrowsWithMessage([]() { parseSpeciesNames(YAML::Load("proton"), species_map); }, "Species is not a sequence!");
  expectThrowsWithMessage(
    []() { parseSpeciesNames(YAML::Load("[proton, electron]"), species_map); },
    "Species was not created in the Species block!");
}

TEST(ParseHelpers, RequiredPositiveDouble) {
  EXPECT_DOUBLE_EQ(parseRequiredPositiveDouble(YAML::Load("Temperature: 2.5"), "Temperature"), 2.5);
  expectThrowsWithMessage(
    []() { parseRequiredPositiveDouble(YAML::Load("Other: 1"), "Temperature"); }, "Temperature is required!");
  expectThrowsWithMessage(
    []() { parseRequiredPositiveDouble(YAML::Load("Temperature: 0."), "Temperature"); }, "Temperature is nonpositive!");
}

TEST(ParseHelpers, VectorMustBeASequence) {
  const mfem::Vector vector = parseVector(YAML::Load("[1., 2.]"), "Center");
  ASSERT_EQ(vector.Size(), 2);
  EXPECT_DOUBLE_EQ(vector[0], 1.);
  EXPECT_DOUBLE_EQ(vector[1], 2.);
  expectThrowsWithMessage([]() { parseVector(YAML::Load("1."), "Center"); }, "Center must be a sequence.");
}

TEST(ParseHelpers, Vector3MustHaveLengthThree) {
  const mfem::Vector vector = parseVector3(YAML::Load("[1., 2., 3.]"), "Bulk Velocity");
  ASSERT_EQ(vector.Size(), 3);
  EXPECT_DOUBLE_EQ(vector[2], 3.);
  const std::string message = "Bulk Velocity does not have length 3!";
  expectThrowsWithMessage([]() { parseVector3(YAML::Load("[1., 2.]"), "Bulk Velocity"); }, message);
  expectThrowsWithMessage([]() { parseVector3(YAML::Load("1."), "Bulk Velocity"); }, message);
}

} // namespace
