#pragma once

#include <libmfpic/DGEulerAssembly.hpp>
#include <libmfpic/Euler.hpp>
#include <mfem.hpp>

namespace mfpic {

/**
 * @brief mfem::LinearFormIntegrator that can interface with \ref ElectromagneticFieldsEvaluator to compute
 * source terms for the Euler-Maxwell system.
 */

class EulerMaxwellSourceIntegrator : public mfem::LinearFormIntegrator 
{

  public:

  /**
   * @brief Construct a new integrator
   *
   * @param fluid_grid_function fluid GridFunction
   * @param field_evaluator electromagnetic field evaluator
   * @param charge_over_mass charge-to-mass ratio for the species
   * @param include_energy_source option to include the total energy source term
   */
  EulerMaxwellSourceIntegrator(const mfem::GridFunction &fluid_grid_function, 
                               const ElectromagneticFieldsEvaluator &field_evaluator, 
                               const mfem::real_t charge_over_mass, 
                               const bool include_energy_source = false) :
    mfem::LinearFormIntegrator(),
    fluid_grid_function_(fluid_grid_function),
    field_evaluator_(field_evaluator),
    charge_over_mass_(charge_over_mass),
    num_equations_(fluid_grid_function.VectorDim()),
    include_energy_source_(include_energy_source) {}

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
    const int element_id = element_transformation.ElementNo;

    rhs.SetSize(num_dof * num_equations_);
    rhs = 0.;
    mfem::DenseMatrix rhs_in_element(rhs.GetData(), num_dof, num_equations_);

    const mfem::IntegrationRule *integration_rule = this->GetIntegrationRule(element, element_transformation);
    if (integration_rule == NULL) {
      integration_rule = &mfem::IntRules.Get(element.GetGeomType(), 2*element.GetOrder() + 1);
    }

    mfem::DenseMatrix fluid_state_at_integration_point_locations, integration_point_locations_in_physical_frame;

    element_transformation.Transform(*integration_rule, integration_point_locations_in_physical_frame);
    fluid_grid_function_.GetVectorValues(element_transformation, *integration_rule, fluid_state_at_integration_point_locations);

    mfem::Vector position(integration_point_locations_in_physical_frame.NumRows());
    mfem::Vector fluid_state(num_equations_);
    mfem::Vector basis_values(num_dof);

    for (int ipoint = 0; ipoint < integration_rule->GetNPoints(); ++ipoint) {

      const mfem::IntegrationPoint &integration_point = integration_rule->IntPoint(ipoint); 
      element_transformation.SetIntPoint(&integration_point);
      element.CalcShape(integration_point, basis_values);

      integration_point_locations_in_physical_frame.GetColumn(ipoint, position);
      fluid_state_at_integration_point_locations.GetColumn(ipoint, fluid_state);

      const mfem::Vector e_field = field_evaluator_.getEFieldAt(position, element_id);
      const mfem::Vector b_field = field_evaluator_.getBFieldAt(position, element_id);
      const mfem::real_t density = fluid_state(euler::ConservativeVariables::MASS_DENSITY);
      const mfem::Vector velocity = euler::getBulkVelocityFromConservativeState(fluid_state);

      const mfem::real_t e_dot_v = e_field * velocity;
      mfem::Vector v_cross_b;
      velocity.cross3D(b_field, v_cross_b);

      const double weight = integration_point.weight * element_transformation.Weight();

      for (int dof = 0; dof < num_dof; ++dof) {

        rhs_in_element(dof, euler::ConservativeVariables::X_MOMENTUM_DENSITY) += weight * basis_values(dof) * 
          charge_over_mass_ * density * (e_field(0) + v_cross_b(0));
        rhs_in_element(dof, euler::ConservativeVariables::Y_MOMENTUM_DENSITY) += weight * basis_values(dof) * 
          charge_over_mass_ * density * (e_field(1) + v_cross_b(1));
        rhs_in_element(dof, euler::ConservativeVariables::Z_MOMENTUM_DENSITY) += weight * basis_values(dof) * 
          charge_over_mass_ * density * (e_field(2) + v_cross_b(2));

        rhs_in_element(dof, euler::ConservativeVariables::TOTAL_ENERGY_DENSITY) += 
          weight * basis_values(dof) * charge_over_mass_ * density * e_dot_v * include_energy_source_;
 
      }
    }
  }

  private:
    const mfem::GridFunction &fluid_grid_function_;
    const ElectromagneticFieldsEvaluator &field_evaluator_;
    const mfem::real_t charge_over_mass_;
    const int num_equations_;
    const bool include_energy_source_;

};

} // namespace
