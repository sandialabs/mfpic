#pragma once

#include <libmfpic/Errors.hpp>
#include <mfem/fem/fespace.hpp>
#include <mfem/mfem.hpp>

namespace mfpic {

/**
 * Supported finite element types
 */
enum class FETypes {
  /// H1 continuous, HGRAD
  HGRAD,
  /// Discontinous-Galerkin, HGRAD (on an element)     
  DG
};

class Discretization {
  public:
  /**
  * @brief Construct a new Discretization object 
  *
  * @param mesh  - the mesh for the discretization
  * @param order - the order of the finite elements
  * @param element_type - see \ref FETypes for support finite element types
  * @param vector_dim - if greater than 1, create a vector finite element space with dimension \p vector_dim
  */
  Discretization(mfem::Mesh *mesh, int order, FETypes element_type = FETypes::HGRAD, int vector_dim = 1)
    : element_type_(element_type)
  {
    switch (element_type) {
      case FETypes::HGRAD:
        finite_element_collection_ = std::make_unique<mfem::H1_FECollection>(order,mesh->Dimension());
        break;
        break;
      case FETypes::DG:
        finite_element_collection_ = std::make_unique<mfem::DG_FECollection>(order,mesh->Dimension());
        break;
      default:
        errorWithDeveloperMessage("Finite element type not supported!");
        break;
    } 
    
    finite_element_space_ = std::make_unique<mfem::FiniteElementSpace>(mesh,finite_element_collection_.get(), vector_dim, mfem::Ordering::byNODES);
  }
  
  /**
  * @brief Get the finite element space associated with discretization
  *
  * @return mfem::FiniteElementSpace - reference to the finite element space 
  */
  mfem::FiniteElementSpace &getFeSpace() {return *finite_element_space_;}
  
  /**
  * @brief Get the finite element space associated with discretization
  *
  * @return mfem::FiniteElementSpace - reference to the finite element space 
  */
  const mfem::FiniteElementSpace &getFeSpace() const {return *finite_element_space_;}
  
  /**
  * @brief Get the finite element collection associated with discretization
  *
  * @return mfem::FiniteElementCollection - reference to the finite element collection 
  */
  mfem::FiniteElementCollection &getFeCollection() {return *finite_element_collection_;}
  
  /**
  * @brief Get the finite element collection associated with discretization
  *
  * @return mfem::FiniteElementCollection - reference to the finite element collection 
  */
  const mfem::FiniteElementCollection &getFeCollection() const {return *finite_element_collection_;}

  FETypes getElementType() const {return element_type_;};

private:
  ///unique pointer to the finite element space on the mesh 
  std::unique_ptr<mfem::FiniteElementSpace> finite_element_space_;
  
  ///unique pointer to the finite element collection on the mesh 
  std::unique_ptr<mfem::FiniteElementCollection> finite_element_collection_;

  ///Element type
  const FETypes element_type_;
};

} // namespace mfpic
