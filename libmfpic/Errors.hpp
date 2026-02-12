#pragma once

#include <source_location>
#include <string>

namespace YAML {
class Node;
}

namespace mfpic {

[[noreturn]] void errorWithDeveloperMessage(
  const std::string message,
  const std::source_location location = std::source_location::current()
);

[[noreturn]] void errorWithUserMessage(const std::string message);

std::string formatParseMessage(const YAML::Node& node, const std::string& message);

} // namespace mfpic
