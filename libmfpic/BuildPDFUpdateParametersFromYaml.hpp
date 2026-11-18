#pragma once

#include <libmfpic/PDFUpdateParameters.hpp>

namespace YAML {
class Node;
}

namespace mfpic {

/// @brief Construct PDFUpdateParameters from YAML.
///
/// Expects an optional top-level "PDF Update" node with an optional
/// "Stride" key (int, default 1).  If the node is null/undefined, returns
/// the default (disabled) parameters.
///
/// @param node YAML node (may be null/undefined).
/// @return PDFUpdateParameters
PDFUpdateParameters buildPDFUpdateParametersFromYaml(const YAML::Node& node);

} // namespace mfpic
