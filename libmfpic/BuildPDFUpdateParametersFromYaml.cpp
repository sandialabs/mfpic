#include <libmfpic/BuildPDFUpdateParametersFromYaml.hpp>
#include <libmfpic/Errors.hpp>

#include <yaml-cpp/yaml.h>

namespace mfpic {

PDFUpdateParameters buildPDFUpdateParametersFromYaml(const YAML::Node& node) {
  PDFUpdateParameters parameters;
  if (!node.IsDefined() || node.IsNull()) return parameters;

  parameters.enabled = true;

  if (node["Stride"]) {
    const int stride = node["Stride"].as<int>();
    if (stride < 1) {
      errorWithUserMessage(
        formatParseMessage(node["Stride"], "PDF Update Stride must be greater than 0!"));
    }
    parameters.stride = stride;
  }

  return parameters;
}

} // namespace mfpic
