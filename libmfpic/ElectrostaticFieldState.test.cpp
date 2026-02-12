#include <libmfpic/ElectrostaticFieldState.hpp>

#include <mfem/fem/fespace.hpp>
#include <mfem/mfem.hpp>

#include <gtest/gtest.h>

namespace {

using namespace mfpic;

TEST(ElectrostaticFieldState, getBFieldAtReturnZero) {
  constexpr int num_elems = 10;
  constexpr mfem::real_t length = 1.0;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems, length);

  constexpr int basis_order = 1;
  Discretization es_discretization(&mesh, basis_order);

  ElectrostaticFieldState es_field_state(es_discretization);

  constexpr double dx = length / num_elems;

  mfem::Vector position{0.5 * dx};
  constexpr int element_index = 0;
  mfem::Vector b_field = es_field_state.getBFieldAt(position, element_index);

  EXPECT_EQ(b_field.Size(), 3);
  EXPECT_EQ(b_field[0], 0.);
  EXPECT_EQ(b_field[1], 0.);
  EXPECT_EQ(b_field[2], 0.);
}

mfem::GridFunction constructLinearGridFunction(
  mfem::FiniteElementSpace& finite_element_space,
  const double slope,
  const double y_intercept)
{
  mfem::GridFunction grid_function(&finite_element_space);

  mfem::FunctionCoefficient line_function([slope, y_intercept](const mfem::Vector& x) -> mfem::real_t {
    return slope * x[0] + y_intercept;
  });

  grid_function.ProjectCoefficient(line_function);

  return grid_function;
}

TEST(ElectrostaticFieldState, getEFieldAtNode) {
  constexpr int num_elems = 10;
  constexpr mfem::real_t length = 1.0;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems, length);

  constexpr int basis_order = 1;
  Discretization es_discretization(&mesh, basis_order);

  ElectrostaticFieldState es_field_state(es_discretization);

  constexpr double slope = 2.3;
  constexpr double y_intercept = 4.5;
  mfem::GridFunction potential = constructLinearGridFunction(es_discretization.getFeSpace(), slope, y_intercept);
  es_field_state.setPotential(potential);

  constexpr int node_index = 1;
  mfem::real_t* vertex_position = mesh.GetVertex(node_index);
  mfem::Vector node_position{vertex_position[0]};

  int inner_element, outer_element;
  mesh.GetFaceElements(node_index, &inner_element, &outer_element);

  mfem::Vector e_field_inner = es_field_state.getEFieldAt(node_position, inner_element);
  mfem::Vector e_field_outer = es_field_state.getEFieldAt(node_position, outer_element);

  EXPECT_EQ(e_field_inner.Size(), 3);
  EXPECT_EQ(e_field_outer.Size(), 3);

  constexpr double absolute_tolerance = 1e-14;
  EXPECT_NEAR(e_field_inner[0], -slope, absolute_tolerance);
  EXPECT_NEAR(e_field_outer[0], -slope, absolute_tolerance);
  EXPECT_EQ(e_field_inner[1], 0.);
  EXPECT_EQ(e_field_outer[1], 0.);
  EXPECT_EQ(e_field_inner[2], 0.);
  EXPECT_EQ(e_field_outer[2], 0.);
}

TEST(ElectrostaticFieldState, getEFieldAtRandomPointInElement) {
  constexpr int num_elems = 10;
  constexpr mfem::real_t length = 1.0;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems, length);

  const int basis_order = 1;
  Discretization es_discretization(&mesh, basis_order);
  ElectrostaticFieldState es_field_state(es_discretization);

  constexpr double slope = 6.7;
  constexpr double y_intercept = 8.9;
  mfem::GridFunction potential = constructLinearGridFunction(es_discretization.getFeSpace(), slope, y_intercept);
  es_field_state.setPotential(potential);

  constexpr double dx = length / num_elems;
  std::mt19937 gen;
  std::uniform_real_distribution<> dis(0., dx);

  constexpr int element_index = 0;
  mfem::Vector position{dis(gen)};
  mfem::Vector e_field = es_field_state.getEFieldAt(position, element_index);

  EXPECT_EQ(e_field.Size(), 3);

  EXPECT_DOUBLE_EQ(e_field[0], -slope);
  EXPECT_EQ(e_field[1], 0.);
  EXPECT_EQ(e_field[2], 0.);
}

}
