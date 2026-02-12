#pragma once

#include <libmfpic/DGEulerAssembly.hpp>
#include <libmfpic/Euler.hpp>
#include <mfem.hpp>

namespace mfpic {

/**
 * @brief mfem::LinearFormIntegrator that computes the integrated kinetic energy \f$\int_{K_i} \frac{1}{2\rho} \bm{p}\cdot\bm{p} \phi \, dV\f$.
 */

class EulerKineticEnergyIntegrator : public mfem::LinearFormIntegrator
{

  public:

  /**
   * @brief Construct a new integrator
   *
   * @param fluid_grid_function fluid GridFunction
   */
  EulerKineticEnergyIntegrator(const mfem::GridFunction &fluid_grid_function) :
    mfem::LinearFormIntegrator(),
    fluid_grid_function_(fluid_grid_function),
    num_equations_(fluid_grid_function.VectorDim()) {}

  /**
   * @brief computes the local rhs contribution on a given element
   *
   * @param element current finite element
   * @param element_transformation current element transformation
   * @param rhs storage for local contribution
   */
  void AssembleRHSElementVect(const mfem::FiniteElement &element, mfem::ElementTransformation &element_transformation, mfem::Vector &rhs) override
  {
    const int num_dof = element.GetDof();

    rhs.SetSize(num_dof * num_equations_);
    rhs = 0.;
    mfem::DenseMatrix rhs_in_element(rhs.GetData(), num_dof, num_equations_);

    const mfem::IntegrationRule *integration_rule = this->GetIntegrationRule(element, element_transformation);
    if (integration_rule == NULL) {
      integration_rule = &mfem::IntRules.Get(element.GetGeomType(), 2*element.GetOrder() + 1);
    }

    mfem::DenseMatrix fluid_state_at_integration_point_locations;

    fluid_grid_function_.GetVectorValues(element_transformation, *integration_rule, fluid_state_at_integration_point_locations);

    mfem::Vector fluid_state(num_equations_);
    mfem::Vector basis_values(num_dof);

    for (int ipoint = 0; ipoint < integration_rule->GetNPoints(); ++ipoint) {

      const mfem::IntegrationPoint &integration_point = integration_rule->IntPoint(ipoint);
      element_transformation.SetIntPoint(&integration_point);
      element.CalcShape(integration_point, basis_values);

      fluid_state_at_integration_point_locations.GetColumn(ipoint, fluid_state);

      const mfem::Vector momentum = euler::getMomentumDensityFromConservativeState(fluid_state);
      const mfem::real_t kinetic_energy = euler::getKineticEnergyDensityFromConservativeState(fluid_state);

      const double weight = integration_point.weight * element_transformation.Weight();

      for (int dof = 0; dof < num_dof; ++dof) {

        rhs_in_element(dof, euler::ConservativeVariables::TOTAL_ENERGY_DENSITY) +=
          weight * basis_values(dof) * kinetic_energy;

      }
    }
  }

  private:
    const mfem::GridFunction &fluid_grid_function_;
    const int num_equations_;

};

} // namespace
