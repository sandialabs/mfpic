#include <libmfpic/MFEMHelpers.hpp>

namespace mfpic {

TransformedVectorCoefficient::TransformedVectorCoefficient(
  std::unique_ptr<mfem::VectorCoefficient>&& vector_coefficient,
  std::function<double(const mfem::Vector&)> transformation)
  : vector_coefficient_(std::move(vector_coefficient))
  , transformation_(transformation)
  {}

void TransformedVectorCoefficient::SetTime(double time) {
  vector_coefficient_->SetTime(time);
  this->Coefficient::SetTime(time);
}

double TransformedVectorCoefficient::Eval(
  mfem::ElementTransformation& element_transformation,
  const mfem::IntegrationPoint& integration_point)
{
  mfem::Vector vec(vector_coefficient_->GetVDim());
  vector_coefficient_->Eval(vec, element_transformation, integration_point);
  return transformation_(vec);
}

}