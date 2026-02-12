#pragma once

#include <mfem/mfem.hpp>

namespace mfpic {

/**
 * @brief Get the centroid of an element's local face.
 *
 * @param[in] mesh         Mesh.
 * @param[in] element      Index of element to query.
 * @param[in] element_face Element-local index of face to query.
 *
 * @returns Centroid of @a element_face on @a element.
 */
mfem::Vector getElementFaceCentroid(const mfem::Mesh& mesh, int element, int element_face);

/**
 * @brief Get the index of the element on the other side of an element's local face.
 *
 * @param[in] mesh         Mesh.
 * @param[in] element      Index of element to query.
 * @param[in] element_face Element-local index of face to query.
 *
 * @returns Index of element on the other side of @a element_face from @a element.
 */
int getElementOnOtherSideOfFace(const mfem::Mesh& mesh, int element, int element_face);

/**
 * @brief Get the unit normal of an element's local face pointing out of the element.
 *
 * @param[in] mesh         Mesh.
 * @param[in] element      Index of element to query.
 * @param[in] element_face Element-local index of face to query.
 *
 * @returns Unit normal of @a element_face pointing away from @a element.
 */
mfem::Vector getElementFaceOutwardUnitNormal(const mfem::Mesh& mesh, int element, int element_face);

/**
 * @brief Get the element and element-local face corresponding to a boundary element (which is a face).
 *
 * Meshes that have been modified to be periodic retain boundary elements that no longer connect to the mesh.
 * This function therefore also returns whether a matching face was found.
 *
 * @param[in] mesh             Mesh.
 * @param[in] boundary_element Index of boundary element.
 *
 * @returns Tuple consisting of (element index, element-local face index, bool of whether the face was found).
 */
std::tuple<int, int, bool> getElementFaceOfBoundaryElement(const mfem::Mesh& mesh, int boundary_element);

/**
 * @brief Compute the volume of the mesh.
 *
 * @todo mfem::Mesh::GetElementVolume is currently non-const. Convince MFEM to add a const version.
 *
 * @param[in] mesh Mesh.
 *
 * @returns Volume of the mesh.
 */
double getMeshVolume(mfem::Mesh& mesh);

/**
 * @brief Get the number of faces on an element of the mesh.
 *
 * @param[in] mesh    Mesh.
 * @param[in] element Element to query.
 *
 * @returns Number of faces on element.
 */
int getNumFacesOnElement(const mfem::Mesh& mesh, int element);

/**
 * @brief Get a map from sides on an inline mesh to the boundary attribute set on those sides
 *  side names are "left", "right" in 1D
 *    "left", "right", "bottom", "top" in 2D
 *    "left", "right", "bottom", "top", "back", "front" in 3D
 *
 * @param mesh_dimension - dimension of the inline mesh
 * @return std::unordered_map<std::string, int> - map from sides names to boundary attributes
 */
std::unordered_map<std::string, int> getSideNameToBoundaryAttributeForInlineMeshes(const int mesh_dimension);

/**
 * @brief Compute the smallest cell lengthscale (e.g., \f$h_{min}\f$) over the mesh.
 *
 * @todo mfem::Mesh::GetElementSize is currently non-const. Convince MFEM to add a const version.
 *
 * @param[in] mesh Mesh.
 *
 * @returns Lengthscale associated with smallest cell
 */
double getSmallestCellLengthscale(mfem::Mesh& mesh);

} // namespace mfpic
