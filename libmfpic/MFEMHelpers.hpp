#pragma once

#include <mfem/mfem.hpp>

namespace mfpic {

class TransformedVectorCoefficient: public mfem::Coefficient {
public:
  TransformedVectorCoefficient(
    std::unique_ptr<mfem::VectorCoefficient>&& vector_coefficient,
    std::function<double(const mfem::Vector&)> transformation);

  virtual void SetTime(double time);

  virtual double Eval(
    mfem::ElementTransformation& element_transformation,
    const mfem::IntegrationPoint& integration_point);

private:
  std::unique_ptr<mfem::VectorCoefficient> vector_coefficient_;
  std::function<double(const mfem::Vector&)> transformation_;
};

}