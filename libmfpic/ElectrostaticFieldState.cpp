#include <libmfpic/ElectrostaticFieldState.hpp>

namespace mfpic {

ElectrostaticFieldState::ElectrostaticFieldState(
  Discretization& electrostatic_discretization)
  : potential_(&electrostatic_discretization.getFeSpace())
{
  potential_ = 0.;
}

mfem::Vector ElectrostaticFieldState::getEFieldAt(const mfem::Vector& position, const int element_index) const {
  mfem::Vector e_field(3);
  e_field = 0.;

  const mfem::FiniteElementSpace* es_finite_element_space = potential_.FESpace();
  const mfem::Mesh* mesh = es_finite_element_space->GetMesh();
  mfem::IsoparametricTransformation element_transformation;
  mesh->GetElementTransformation(element_index, &element_transformation);

  mfem::IntegrationPoint integration_point;
  element_transformation.TransformBack(position, integration_point);
  element_transformation.SetIntPoint(&integration_point);

  mfem::Vector gradient;
  potential_.GetGradient(element_transformation, gradient);

  for (int i = 0; i < gradient.Size(); ++i) {
    e_field[i] = -gradient[i];
  }

  return e_field;
}

mfem::Vector ElectrostaticFieldState::getBFieldAt(const mfem::Vector&, const int) const {
  mfem::Vector b_field(3);
  b_field = 0.;
  return b_field;
}



}