#pragma once

#include <libmfpic/DGBC.hpp>
#include <libmfpic/DGGhostBoundaryIntegrator.hpp>
#include <mfem.hpp>

namespace mfpic {

/**
 * @brief A DG BC that uses the ghost cell approach.
 */
struct DGGhostBC : public DGBC {

  DGGhostBC() = delete;
  DGGhostBC(const int boundary_attribute, const mfem::Mesh& mesh, std::unique_ptr<GhostDOFSetter> && dof_setter) 
  : DGBC(boundary_attribute, mesh), ghost_dof_setter(std::move(dof_setter)) {};

  virtual ~DGGhostBC() = default;

  std::shared_ptr<mfem::HyperbolicFormIntegrator> makeIntegrator(const mfem::NumericalFlux & numerical_flux, Species) override {
    return std::make_shared<DGGhostBoundaryIntegrator>(numerical_flux, *ghost_dof_setter);
  }

  std::unique_ptr<GhostDOFSetter> ghost_dof_setter;

};

} // namespace
