#include <libmfpic/Errors.hpp>

#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include <stdexcept>

namespace {

using namespace mfpic;

TEST(Errors, DeveloperMessages) {
  std::string message = "test";

  ASSERT_DEATH(errorWithDeveloperMessage(message), message);
}

TEST(Errors, UserMessages) {
  std::string message = "test";

  ASSERT_THROW(errorWithUserMessage(message), std::logic_error);
}

TEST(Errors, FormatParseMessage) {
  const std::string yaml("Line One Column Zero Node:");
  YAML::Node main = YAML::Load(yaml);

  const std::string inserted_message = "inserted";
  const std::string formatted_parse_message = formatParseMessage(main, inserted_message);

  EXPECT_NE(formatted_parse_message.find(inserted_message), std::string::npos);
  EXPECT_NE(formatted_parse_message.find("line 1, column 0"), std::string::npos);
}

} // namespace
