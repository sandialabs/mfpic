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

std::pair<mfem::Vector, int> MeshDistribution::generateRandomPointAndElement(RandomNumberGenerator& generator) {
  const int random_element = element_distribution_(generator);

  mfem::IntegrationPoint random_reference_point;
  const mfem::Geometry::Type element_geometry = mesh_->GetElementGeometry(random_element);
  // TODO MFEM currently calls rand() to generate random numbers. oughtta use c++ generators.
  // The following is done to provide a seed consistent with our own.
  std::uniform_int_distribution<unsigned> random_unsigned_distribution;
  std::srand(random_unsigned_distribution(generator));
  mfem::Geometry::GetRandomPoint(element_geometry, random_reference_point);

  mfem::Vector random_physical_point;
  mfem::ElementTransformation* element_transformation = mesh_->GetElementTransformation(random_element);
  element_transformation->Transform(random_reference_point, random_physical_point);

  return std::make_pair(random_physical_point, random_element);
}

} // namespace mfpic
