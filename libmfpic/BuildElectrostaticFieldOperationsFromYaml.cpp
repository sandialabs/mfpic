#include <libmfpic/BuildElectrostaticFieldOperationsFromYaml.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/ElectrostaticFieldOperationsWithBoltzmannElectrons.hpp>

#include <yaml-cpp/yaml.h>

namespace mfpic {

std::unique_ptr<ElectrostaticFieldOperations> buildElectrostaticFieldOperations(
  const YAML::Node& fields_node,
  Discretization& electrostatic_discretization,
  std::unique_ptr<DirichletBoundaryConditions> dirichlet_boundary_conditions
) {
  std::unique_ptr<ElectrostaticFieldOperations> electrostatic_field_operations;

  if (fields_node["Boltzmann Electrons"]) {
    const YAML::Node& boltzmann_electrons_node = fields_node["Boltzmann Electrons"];
    const double reference_number_density = boltzmann_electrons_node["Reference Number Density"].as<double>();
    const double temperature = boltzmann_electrons_node["Temperature"].as<double>();

    electrostatic_field_operations = std::make_unique<ElectrostaticFieldOperationsWithBoltzmannElectrons>(
      electrostatic_discretization,
      std::move(dirichlet_boundary_conditions),
      reference_number_density,
      temperature
    );
  } else {
    electrostatic_field_operations = std::make_unique<ElectrostaticFieldOperations>(
      electrostatic_discretization,
      std::move(dirichlet_boundary_conditions)
    );
  }

  return electrostatic_field_operations;
}

} // namespace mfpic
