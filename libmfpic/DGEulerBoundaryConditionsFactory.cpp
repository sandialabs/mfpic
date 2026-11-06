#include <libmfpic/DGEulerBoundaryConditions.hpp>
#include <libmfpic/DGEulerBoundaryConditionsFactory.hpp>
#include <libmfpic/DGGhostBC.hpp>
#include <libmfpic/Errors.hpp>
#include <libmfpic/KineticFluxBC.hpp>
#include <libmfpic/MaxwellianInfluxBC.hpp>
#include <libmfpic/ParseHelpers.hpp>

#include <yaml-cpp/yaml.h>

#include <functional>
#include <memory>

namespace mfpic {

std::unique_ptr<DGBC> ReflectingDGEulerBCParameters::makeBC(
  int boundary_attribute,
  const mfem::Mesh& mesh,
  const Species& species
) const {
  auto dof_setter = std::make_unique<DGEulerReflectingBC>();
  return std::make_unique<DGGhostBC>(boundary_attribute, mesh, species, std::move(dof_setter));
}

std::unique_ptr<DGBC> AbsorbingDGEulerBCParameters::makeBC(
  int boundary_attribute,
  const mfem::Mesh& mesh,
  const Species& species
) const {
  return std::make_unique<KineticFluxBC>(boundary_attribute, mesh, species);
}

MaxwellianInfluxDGEulerBCParameters::MaxwellianInfluxDGEulerBCParameters(double number_flux_density_in, double temperature_in)
  : number_flux_density(number_flux_density_in), temperature(temperature_in)
{}

std::unique_ptr<DGBC> MaxwellianInfluxDGEulerBCParameters::makeBC(
  int boundary_attribute,
  const mfem::Mesh& mesh,
  const Species& species
) const {
  return std::make_unique<MaxwellianInfluxBC>(boundary_attribute, mesh, species, number_flux_density, temperature);
}

namespace {

using DGEulerBCParametersParser = std::function<std::unique_ptr<DGEulerBCParameters>(const YAML::Node&)>;

std::unique_ptr<DGEulerBCParameters> parseMaxwellianInfluxBCParameters(const YAML::Node& fluid_bc)
{
  const double number_flux_density = parseRequiredPositiveDouble(fluid_bc, "Number Flux Density");
  const double temperature = parseRequiredPositiveDouble(fluid_bc, "Temperature");
  return std::make_unique<MaxwellianInfluxDGEulerBCParameters>(number_flux_density, temperature);
}

std::unordered_map<std::string, DGEulerBCParametersParser> getBCTypeNameToParser()
{
  std::unordered_map<std::string, DGEulerBCParametersParser> bc_type_name_to_parser{
    {"Reflecting", [](const YAML::Node&) { return std::make_unique<ReflectingDGEulerBCParameters>(); }},
    {"Absorbing", [](const YAML::Node&) { return std::make_unique<AbsorbingDGEulerBCParameters>(); }},
    {"Maxwellian Influx", parseMaxwellianInfluxBCParameters},
  };
  return bc_type_name_to_parser;
}

}

std::unordered_map<int, std::unique_ptr<DGEulerBCParameters>> buildBoundaryAttributeToBCParametersFromYAML(
  const YAML::Node& fluid_bcs,
  const int mesh_dimension)
{
  assert(fluid_bcs.IsSequence() or fluid_bcs.IsNull() or not fluid_bcs.IsDefined());

  std::unordered_map<int, std::unique_ptr<DGEulerBCParameters>> boundary_attribute_to_bc_parameters;
  for (const YAML::Node& fluid_bc : fluid_bcs) {
    const int boundary_attribute = parseBoundaryAttribute(fluid_bc, mesh_dimension);

    const std::string bc_type_string = fluid_bc["Type"].as<std::string>();
    const auto& bc_type_name_to_parser = getBCTypeNameToParser();
    if (not bc_type_name_to_parser.contains(bc_type_string)) {
      errorWithUserMessage(formatParseMessage(fluid_bc["Type"], "This boundary condition type is invalid."));
    }
    boundary_attribute_to_bc_parameters[boundary_attribute] = bc_type_name_to_parser.at(bc_type_string)(fluid_bc);
  }
  return boundary_attribute_to_bc_parameters;
}

std::vector<std::vector<std::unique_ptr<DGBC>>> buildDGEulerBoundaryConditions(
  const std::unordered_map<int, std::unique_ptr<DGEulerBCParameters>>& boundary_attribute_to_bc_parameters,
  const mfem::Mesh& mesh,
  const std::vector<Species>& species_list)
{
  std::vector<std::vector<std::unique_ptr<DGBC>>> dg_euler_bcs;

  for (const Species& species : species_list) {
    std::vector<std::unique_ptr<DGBC>> species_bcs;
    for (const auto& [boundary_attribute, bc_parameters] : boundary_attribute_to_bc_parameters) {
      species_bcs.push_back(bc_parameters->makeBC(boundary_attribute, mesh, species));
    }
    dg_euler_bcs.push_back(std::move(species_bcs));
  }

  return dg_euler_bcs;
}

}
