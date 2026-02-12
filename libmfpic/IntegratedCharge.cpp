#include <libmfpic/IntegratedCharge.hpp>

namespace mfpic {

IntegratedCharge::IntegratedCharge(Discretization &discretization)
: integrated_charge_(discretization.getFeSpace().GetNDofs())
{
  integrated_charge_ = 0.0;
}

void IntegratedCharge::addCharge(const IntegratedCharge& charge_to_add) {
  integrated_charge_ += charge_to_add.integrated_charge_;
}

void IntegratedCharge::setIntegratedCharge(const mfem::Vector& integrated_charge) {
  integrated_charge_ = integrated_charge;
}

mfem::Vector IntegratedCharge::getIntegratedCharge() const { return integrated_charge_; }

void IntegratedCharge::setIntegratedChargeValue(int dof_id, double value) {
  integrated_charge_(dof_id) = value;
}

void IntegratedCharge::addIntegratedChargeValue(int dof_id, double value) {
  integrated_charge_(dof_id) += value;
}

void IntegratedCharge::setIntegratedChargeValue(double value) {
  integrated_charge_ = value;
}

double IntegratedCharge::getIntegratedChargeValue(int dof_id) {
  return integrated_charge_(dof_id);
}

} // namespace mfpic
