#pragma once

#include <libmfpic/RandomNumberGenerator.hpp>

#include <mfem/mfem.hpp>

namespace mfpic {

/// Generate points in a mesh according to a given probability density function.
class MeshDistribution {
public:
  /**
   * @brief Ctor.
   *
   * @param[in] mesh         Mesh in which to generate points.
   * @param[in] distribution Probability density function of points in the mesh (no need to normalize).
   */
  MeshDistribution(std::shared_ptr<mfem::Mesh> mesh, std::function<double(const mfem::Vector&)> distribution);

  /**
   * @brief Generate a random point in the mesh and the index of the element containing it.
   *
   * @param[in,out] generator Random number generator.
   *
   * @returns Both a point in the mesh having length equal to the mesh dimension and the element index containing it.
   */
  std::pair<mfem::Vector, int> generateRandomPointAndElement(RandomNumberGenerator& generator);

private:
  /// Mesh in which to generate points.
  std::shared_ptr<mfem::Mesh> mesh_;

  /// Probability distribution for the elements, weighting each by its volume.
  std::discrete_distribution<> element_distribution_;
};

} // namespace mfpic
