#include <libmfpic/ElectrostaticFieldState.hpp>

namespace mfpic {

namespace {
constexpr int e_field_dg_order = 0;
constexpr int e_field_dim = 3;
}

ElectrostaticFieldState::ElectrostaticFieldState(
  Discretization& electrostatic_discretization)
  : potential_(&electrostatic_discretization.getFeSpace())
  , e_field_dg_discretization_(
      std::make_shared<Discretization>(
        electrostatic_discretization.getFeSpace().GetMesh(),
        e_field_dg_order,
        FETypes::DG,
        e_field_dim))
  , e_field_(&e_field_dg_discretization_->getFeSpace())
{
  potential_ = 0.;
  e_field_ = 0.;
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