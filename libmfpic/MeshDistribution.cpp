#include <libmfpic/MeshDistribution.hpp>

namespace mfpic {

MeshDistribution::MeshDistribution(std::shared_ptr<mfem::Mesh> mesh) : mesh_(mesh) {
  const int num_elements = mesh->GetNE();
  std::vector<double> element_volumes(num_elements);
  for (int element = 0; element < num_elements; element++) {
    element_volumes[element] = mesh->GetElementVolume(element);
  }

  element_distribution_ = std::discrete_distribution<>(element_volumes.cbegin(), element_volumes.cend());
}

} // namespace mfpic
