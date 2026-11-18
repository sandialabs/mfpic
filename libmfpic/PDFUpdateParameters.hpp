#pragma once

namespace mfpic {

/// @brief Holds options for periodic re-estimation of per-particle f_HF values.
struct PDFUpdateParameters {
  bool enabled = false;
  int stride = 1;
};

} // namespace mfpic
