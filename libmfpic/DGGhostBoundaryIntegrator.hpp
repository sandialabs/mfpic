#pragma once

#include <libmfpic/DGGhostBC.hpp>
#include <mfem/fem/hyperbolic.hpp>
#include <mfem/linalg/densemat.hpp>

namespace mfpic {

/**
 * @brief mfem::HyperbolicFormIntegrator that applies the ghost cell approach.
 */

 class DGGhostBoundaryIntegrator : public mfem::HyperbolicFormIntegrator {

  public:

  /**
   * @brief Construct a new DGGhostBoundaryIntegrator
   * 
   * @param numerical_flux mfem::NumericalFlux function
   * @param ghost_dof_setter a \ref DGGhostDOFSetter for the ghost cell computation 
   */

  DGGhostBoundaryIntegrator(const mfem::NumericalFlux & numerical_flux,
                            const DGGhostBC & ghost_dof_setter) :
    mfem::HyperbolicFormIntegrator(numerical_flux),
    num_equations_(numerical_flux.GetFluxFunction().num_equations),
    ghost_dof_setter_(ghost_dof_setter) {};

  /**
   * @brief A mfem::FaceElementTransformations object for ghost cells. 
   * Copies the 1 data into the 2 data.
   *
   * @note Intended only to be used through copying the transformations
   * from the real, interior boundary-touching cell.
   */
  
  class GhostFaceElementTransformations : public mfem::FaceElementTransformations {

    public:

    GhostFaceElementTransformations() = delete;
    GhostFaceElementTransformations & operator=(const GhostFaceElementTransformations & existing) = delete;
    GhostFaceElementTransformations(const mfem::FaceElementTransformations & transformations) :
      mfem::FaceElementTransformations(transformations) {
        Elem2No = transformations.Elem1No;
        Elem2   = transformations.Elem1;
        Loc2    = transformations.Loc1;
        updateMaskAfterGhostIsSet_();
      };
    ~GhostFaceElementTransformations() = default;

    private:
      void updateMaskAfterGhostIsSet_() {SetConfigurationMask(mask_);};
      const int mask_ = (0 | 1 | 2 | 4 | 8 | 16);
  };

  /**
   * @brief Compute the hyperbolic flux terms on the domain boundary, \f$\left\langle\hat{F} \cdot n, v\right\rangle\f$ 
   * and stores in \p rhs .
   *
   * @param element_left  mfem::FiniteElement for the real, interior cell
   * @param transformations face element transformations
   * @param element_dofs mfem::Vector of element DOFs for the relevant cells
   * @param[out] rhs storage for the integrated flux
   */

  void AssembleFaceVector(const mfem::FiniteElement & element_left,
                          const mfem::FiniteElement &,
                          mfem::FaceElementTransformations & transformations,
                          const mfem::Vector & element_dofs,
                          mfem::Vector & rhs) override 
  {

    const int num_dof = element_left.GetDof(); 

    mfem::Vector element_dofs_with_ghost(2 * num_dof * num_equations_);
    assert(element_dofs.Size() == num_dof * num_equations_);
    element_dofs_with_ghost.SetVector(element_dofs,0);

    const mfem::DenseMatrix interior_dofs_as_matrix(element_dofs_with_ghost.GetData(), num_dof, num_equations_);
    mfem::DenseMatrix ghost_dofs_as_matrix(element_dofs_with_ghost.GetData() + num_dof * num_equations_, num_dof, num_equations_);

    transformations.SetIntPoint(&mfem::Geometries.GetCenter(transformations.GetGeometryType()));

    mfem::Vector normal{0.,0.,0.};
    if (transformations.GetSpaceDim() == 1)
    {
      normal(0) = 2*transformations.GetElement1IntPoint().x - 1.;
    }
    else
    {
      mfem::CalcOrtho(transformations.Jacobian(), normal);
    }
    normal /= normal.Norml2();

    ghost_dof_setter_.setDOFsInGhost(interior_dofs_as_matrix, normal, ghost_dofs_as_matrix);

    GhostFaceElementTransformations ghost_transformations = transformations;

    mfem::HyperbolicFormIntegrator::AssembleFaceVector(element_left, 
                                                       element_left,
                                                       ghost_transformations,
                                                       element_dofs_with_ghost,
                                                       rhs);
    rhs.SetSize(num_dof * num_equations_);
  }

  private:
    /// number of equations in the system
    const int num_equations_;
    const DGGhostBC & ghost_dof_setter_;

};

} // namespace
