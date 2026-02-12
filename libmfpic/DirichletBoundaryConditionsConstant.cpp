#include <libmfpic/DirichletBoundaryConditionsConstant.hpp>

namespace mfpic {

DirichletBoundaryConditionsConstant::DirichletBoundaryConditionsConstant(
  const std::unordered_map<int, double>& boundary_attribute_to_dirichlet_value,
  Discretization& discretization)
  : boundary_attribute_to_dirichlet_value_(boundary_attribute_to_dirichlet_value)
{
  mfem::FiniteElementSpace& finite_element_space = discretization.getFeSpace();

  const mfem::Mesh* mesh = finite_element_space.GetMesh();
  mfem::Array<int> boundary_attribute_is_dirichlet(mesh->bdr_attributes.Max());
  boundary_attribute_is_dirichlet = false;
  for (const auto& [boundary_attribute, _] : boundary_attribute_to_dirichlet_value_) {
    boundary_attribute_is_dirichlet[boundary_attribute - 1] = true;
  }

  finite_element_space.GetEssentialTrueDofs(boundary_attribute_is_dirichlet, dirichlet_dof_indices_);

}

mfem::Array<int> DirichletBoundaryConditionsConstant::getDirichletBoundaryDofIndices() const { return dirichlet_dof_indices_; }

void DirichletBoundaryConditionsConstant::applyBoundaryConditions(mfem::GridFunction& solution) const {
  for (const auto & [boundary_attribute, dirichlet_value] : boundary_attribute_to_dirichlet_value_) {
    mfem::ConstantCoefficient constant_coefficient(dirichlet_value);

    const mfem::Mesh* mesh = solution.FESpace()->GetMesh();
    mfem::Array<int> project_onto_boundary_attribute(mesh->bdr_attributes.Max());
    project_onto_boundary_attribute = false;
    project_onto_boundary_attribute[boundary_attribute - 1] = true;

    solution.ProjectBdrCoefficient(constant_coefficient, project_onto_boundary_attribute);
  }
}

}
