#pragma once

#include <mfem/mfem.hpp>

#include <random>

namespace mfpic {

/// Generate points uniformly throughout a mesh.
class UniformMeshDistribution {
public:
  /**
   * @brief Ctor.
   *
   * @param[in] mesh Mesh in which to generate points.
   */
  UniformMeshDistribution(std::shared_ptr<mfem::Mesh> mesh);

  /**
   * @brief Generate a random point in the mesh and the index of the element containing it.
   *
   * @tparam Generator A UniformRandomBitGenerator type.
   *
   * @param[in,out] generator A UniformRandomBitGenerator used to generate some random numbers.
   *
   * @returns Both a point in the mesh having length equal to the mesh dimension and the element index containing it.
   */
  template <std::uniform_random_bit_generator Generator>
  std::pair<mfem::Vector, int> generateRandomPointAndElement(Generator& generator);

private:
  /// Mesh in which to generate points.
  std::shared_ptr<mfem::Mesh> mesh_;

  /// Probability distribution for the elements, weighting each by its volume.
  std::discrete_distribution<> element_distribution_;
};

// template definitions
template <std::uniform_random_bit_generator Generator>
std::pair<mfem::Vector, int> UniformMeshDistribution::generateRandomPointAndElement(Generator& generator) {
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
