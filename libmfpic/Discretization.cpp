#include <libmfpic/Discretization.hpp>

namespace mfpic {

Discretization getIdentityTestFunctionDiscretization(mfem::Mesh& mesh) {
  constexpr int finite_element_order = 0;
  constexpr FETypes element_type = FETypes::DG;
  Discretization identity_test_function_discretization(&mesh, finite_element_order, element_type);
  return identity_test_function_discretization;
}

Discretization getIdentityTestFunctionDiscretization(Discretization& discretization){
  mfem::FiniteElementSpace& finite_element_space = discretization.getFeSpace();
  mfem::Mesh* mesh = finite_element_space.GetMesh();
  return getIdentityTestFunctionDiscretization(*mesh);
}

}