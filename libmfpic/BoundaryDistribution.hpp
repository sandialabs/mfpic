#pragma once

#include <libmfpic/RandomNumberGenerator.hpp>

#include <mfem/mfem.hpp>

namespace mfpic {

/// A point on the mesh boundary along with the element it borders and the direction into the domain.
struct BoundaryPoint {
  /// Position of the point, of length 3.
  mfem::Vector position;

  /// Index of the element adjacent to the boundary at the point.
  int element;

  /// Unit normal of the boundary pointing out of the domain, of length 3.
  mfem::Vector outward_unit_normal;
};

/// Uniformly generate points on a boundary of a mesh.
class BoundaryDistribution {
public:
  /**
   * @brief Ctor.
   *
   * @throws If mesh contains no boundary faces with @a boundary_attribute.
   *
   * @param[in] mesh               Mesh whose boundary to sample.
   * @param[in] boundary_attribute Boundary attribute of the faces to sample.
   */
  BoundaryDistribution(std::shared_ptr<mfem::Mesh> mesh, int boundary_attribute);

  /**
   * @brief Get the total measure of the sampled boundary faces.
   *
   * @returns Total area in 3D, length in 2D, and number of boundary points in 1D.
   */
  double getTotalArea() const;

  /**
   * @brief Generate a random point on the indicated boundary of the mesh.
   *
   * @param[in,out] generator Random number generator.
   *
   * @returns A point on the boundary.
   */
  BoundaryPoint generateRandomBoundaryPoint(RandomNumberGenerator& generator);

private:
  /// Data for a face on the boundary.
  struct FaceInfo_ {
    int boundary_element;             ///< Mesh boundary element corresponding to this face.
    int element;                      ///< Mesh element adjacent to this face.
    int element_face;                 ///< Mesh element-local face index corresponding to this face.
    mfem::Vector outward_unit_normal; ///< Outward-facing unit normal.
  };

  /// Mesh whose boundary to sample.
  std::shared_ptr<mfem::Mesh> mesh_;

  /// Data for faces on this boundary.
  std::vector<FaceInfo_> faces_info_;

  /// Total area of this boundary.
  double total_area_ = 0.0;

  /// Probability distribution for the faces, weighting each by its area.
  std::discrete_distribution<> face_distribution_;
};

} // namespace mfpic
