#include <libmfpic/MeshDistribution.hpp>
#include <libmfpic/MeshUtilities.hpp>

namespace mfpic {

MeshDistribution::MeshDistribution(
  std::shared_ptr<mfem::Mesh> mesh,
  std::function<double(const mfem::Vector&)> distribution
) :
  mesh_(mesh)
{
  const mfem::Vector elementwise_distribution_integrals = elementwiseIntegral(*mesh, distribution);
  element_distribution_ = std::discrete_distribution<>(
    elementwise_distribution_integrals.begin(),
    elementwise_distribution_integrals.end()
  );
}

} // namespace mfpic
