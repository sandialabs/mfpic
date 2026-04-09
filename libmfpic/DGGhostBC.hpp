#pragma once

#include <libmfpic/DGBC.hpp>
#include <libmfpic/DGGhostBoundaryIntegrator.hpp>
#include <memory>
#include <mfem.hpp>
#include <mfem/fem/nonlininteg.hpp>

namespace mfpic {

/**
 * @brief A DG BC that uses the ghost cell approach.
 */
struct DGGhostBC : public DGBC {

  DGGhostBC() = delete;
  DGGhostBC(const int boundary_attribute, const mfem::Mesh& mesh, const Species species, std::unique_ptr<GhostDOFSetter> && dof_setter) 
  : DGBC(boundary_attribute, mesh, species), ghost_dof_setter(std::move(dof_setter)) {};

  virtual ~DGGhostBC() = default;

  std::unique_ptr<mfem::NonlinearFormIntegrator> makeIntegrator(const mfem::NumericalFlux & numerical_flux) override {
    return std::make_unique<DGGhostBoundaryIntegrator>(numerical_flux, *ghost_dof_setter);
  }

  std::unique_ptr<GhostDOFSetter> ghost_dof_setter;

};

} // namespace
