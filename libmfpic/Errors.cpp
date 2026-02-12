#include <libmfpic/Errors.hpp>

#include <yaml-cpp/yaml.h>

#include <iostream>

namespace mfpic {

[[noreturn]] void errorWithDeveloperMessage(
  const std::string message,
  const std::source_location location
) {
  std::cerr
  << "file: "
  << location.file_name() << "("
  << location.line() << ":"
  << location.column() << ") `"
  << location.function_name() << "`: "
  << message
  << std::endl;

  std::abort();
}

[[noreturn]] void errorWithUserMessage(const std::string message) {
  throw std::logic_error(message);
}

std::string formatParseMessage(const YAML::Node& node, const std::string& message) {
  const YAML::Mark& mark = node.Mark();
  std::stringstream formatted_message;
  formatted_message << "In input deck, at line " << mark.line + 1 << ", column " << mark.column << ":\n  " << message;
  return formatted_message.str();
}

} // namespace mfpic
