#include <libmfpic/Constants.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/SourcesFactory.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace {

using namespace mfpic;
using namespace testing;

const std::string electron_species_name = "electron";
const std::string proton_species_name = "proton";

const Species electron_species{.charge = -constants::elementary_charge, .mass = constants::electron_mass};
const Species proton_species{.charge = constants::elementary_charge, .mass = constants::proton_mass};

const std::unordered_map<std::string, Species> species_map{
  {electron_species_name, electron_species},
  {proton_species_name, proton_species}};

mfem::Vector evaluateVectorCoefficientAtPoint(
  mfem::VectorCoefficient& coefficient,
  const double x,
  const double dx,
  const mfem::Mesh& mesh)
{
  const int element_index = x / dx;

  mfem::IsoparametricTransformation element_transformation;
  mesh.GetElementTransformation(element_index, &element_transformation);

  const mfem::Vector x_vec{x};
  mfem::IntegrationPoint integration_point;
  element_transformation.TransformBack(x_vec, integration_point);
  element_transformation.SetIntPoint(&integration_point);

  mfem::Vector state_out(euler::ConservativeVariables::NUM_VARS);
  coefficient.Eval(state_out, element_transformation, integration_point);
  return state_out;
}

TEST(SourcesFactory, SodSourceParametersEulerVectorCoefficient) {
  constexpr double discontinuity_location = 0.5;

  SourceStateParameters left_state{.number_density = 1e21, .temperature = 300};
  SourceStateParameters right_state{.number_density = 1e22, .temperature = 320};
  SodSourceParameters parameters(electron_species, discontinuity_location, left_state, right_state);

  std::unique_ptr<mfem::VectorCoefficient> euler_coefficient = parameters.getEulerVectorCoefficient();

  constexpr int num_elems = 10;
  constexpr double length = 1.0;
  constexpr double dx = length / num_elems;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(10);

  const double x_left = parameters.discontinuity_location - 0.5 * dx;
  const mfem::Vector left_state_out = evaluateVectorCoefficientAtPoint(*euler_coefficient, x_left, dx, mesh);

  const mfem::Vector left_state_primitive = euler::constructPrimitiveState(
    left_state.number_density, left_state.bulk_velocity, left_state.temperature);
  const mfem::Vector left_state_expected = euler::convertFromPrimitiveToConservative(left_state_primitive, parameters.species);

  for (int i = 0; i < left_state_out.Size(); ++i) {
    EXPECT_DOUBLE_EQ(left_state_expected[i], left_state_out[i]);
  }

  const double x_right = parameters.discontinuity_location + 0.5 * dx;
  const mfem::Vector right_state_out = evaluateVectorCoefficientAtPoint(*euler_coefficient, x_right, dx, mesh);

  const mfem::Vector right_state_primitive = euler::constructPrimitiveState(
    right_state.number_density, right_state.bulk_velocity, right_state.temperature);
  const mfem::Vector right_state_expected = euler::convertFromPrimitiveToConservative(right_state_primitive, parameters.species);

  for (int i = 0; i < right_state_out.Size(); ++i) {
    EXPECT_DOUBLE_EQ(right_state_expected[i], right_state_out[i]);
  }
}

mfem::Vector evaluateGaussianAtPoint(
  const double x,
  const mfem::Vector& center,
  const double standard_deviation,
  const SourceStateParameters& offsets,
  const SourceStateParameters& heights,
  const double pressure_offset,
  const double pressure_height,
  const Species& species)
{
  const double shift = x - center[0];
  const double exponential = exp(-0.5 * shift * shift / (standard_deviation * standard_deviation));
  const double expected_number_density = heights.number_density * exponential + offsets.number_density;

  mfem::Vector expected_velocity(offsets.bulk_velocity);
  expected_velocity.Add(exponential, heights.bulk_velocity);

  const double expected_pressure = pressure_height * exponential + pressure_offset;

  const double expected_temperature = euler::temperature(expected_number_density, expected_pressure);

  const mfem::Vector expected_state_primitive = euler::constructPrimitiveState(
    expected_number_density, expected_velocity, expected_temperature);

  const mfem::Vector expected_state_conservative = euler::convertFromPrimitiveToConservative(
    expected_state_primitive, species);

  return expected_state_conservative;
}

TEST(SourcesFactory, GaussianSourceParametersEulerVectorCoefficient) {
  const mfem::Vector center{0.5};
  constexpr double standard_deviation = 0.01;

  const mfem::Vector bulk_velocity_offset{1., 0., 0.};
  const mfem::Vector bulk_velocity_height{0.1, 0., 0.};
  SourceStateParameters offsets{.number_density = 1e19, .bulk_velocity=bulk_velocity_offset};
  SourceStateParameters heights{.number_density = 1e18, .bulk_velocity=bulk_velocity_height};

  const double pressure_offset = 0.05;
  const double pressure_height = 0.001;

  GaussianSourceParameters parameters(
    electron_species, center, standard_deviation, offsets, heights, pressure_offset, pressure_height);
  std::unique_ptr<mfem::VectorCoefficient> euler_coefficient = parameters.getEulerVectorCoefficient();

  constexpr int num_elems = 11;
  constexpr double length = 1.0;
  constexpr double dx = length / num_elems;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(10);

  for (const double& x : {0.01, 0.45, 0.5}){
    const mfem::Vector state_out = evaluateVectorCoefficientAtPoint(*euler_coefficient, x, dx, mesh);
    const mfem::Vector expected_state = evaluateGaussianAtPoint(
      x, center, standard_deviation, offsets, heights, pressure_offset, pressure_height, electron_species);
    for (int i = 0; i < expected_state.Size(); ++i) {
      EXPECT_DOUBLE_EQ(expected_state[i], state_out[i]);
    }
  }
}

TEST(SourcesFactory, EmptyMainGivesEmptyListOfSourceParameters) {
  const std::string main_string("");

  YAML::Node main_node = YAML::Load(main_string);
  YAML::Node model_node = main_node["Euler Fluids"];
  YAML::Node initial_conditions_node = model_node["Initial Conditions"];

  const std::vector<std::unique_ptr<SourceParameters>> list_of_parameters = buildListOfSourceParametersFromYAML(
    initial_conditions_node, species_map);

  EXPECT_EQ(0, list_of_parameters.size());
}

TEST(SourcesFactory, NoInitialConditionsGivesEmptyListOfSourceParameters) {
  const std::string main_string("Euler Fluids:");

  YAML::Node main_node = YAML::Load(main_string);
  YAML::Node model_node = main_node["Euler Fluids"];
  YAML::Node initial_conditions_node = model_node["Initial Conditions"];

  const std::vector<std::unique_ptr<SourceParameters>> list_of_parameters = buildListOfSourceParametersFromYAML(
    initial_conditions_node, species_map);

  EXPECT_EQ(0, list_of_parameters.size());
}

TEST(SourcesFactory, EmptyInitialConditionsGivesEmptyListOfSourceParameters) {
  const std::string model_string("Initial Conditions:");

  const YAML::Node model_node = YAML::Load(model_string);
  const YAML::Node initial_conditions_node = model_node["Initial Conditions"];

  const std::vector<std::unique_ptr<SourceParameters>> list_of_parameters = buildListOfSourceParametersFromYAML(
    initial_conditions_node, species_map);

  EXPECT_EQ(0, list_of_parameters.size());
}

TEST(SourcesFactory, SpeciesMustBeASequence) {
  const std::string model_string(
    "Initial Conditions:\n"
    "  - Species: electron\n"
    "    Constant: \n"
    "      Number Density: 12\n"
    "      Temperature: 300\n");

  const YAML::Node model_node = YAML::Load(model_string);
  const YAML::Node initial_conditions_node = model_node["Initial Conditions"];

  EXPECT_THAT(
    [&]() {buildListOfSourceParametersFromYAML(initial_conditions_node, species_map);},
    ThrowsMessage<std::logic_error>(HasSubstr("Species is not a sequence!"))
  );
}

TEST(SourcesFactory, SpeciesMustExistInSpeciesMap) {
  const std::string model_string(
    "Initial Conditions:\n"
    "  - Species: [nonexistent]\n"
    "    Constant:\n"
    "      Number Density: 12\n"
    "      Temperature: 300\n");

  const YAML::Node model_node = YAML::Load(model_string);
  const YAML::Node initial_conditions_node = model_node["Initial Conditions"];

  EXPECT_THAT(
    [&]() {buildListOfSourceParametersFromYAML(initial_conditions_node, species_map);},
    ThrowsMessage<std::logic_error>(HasSubstr("Species was not created in the Species block!"))
  );
}

TEST(SourcesFactory, NumberDensityMustBePositive) {
  const std::string model_string(
    "Constant:\n"
    "  Number Density: 0.\n"
    "  Temperature: 300\n");

  const YAML::Node model_node = YAML::Load(model_string);
  const YAML::Node state_node = model_node["Constant"];

  EXPECT_THAT(
    [&]() {buildSourceStateParametersFromYAML(state_node);},
    ThrowsMessage<std::logic_error>(HasSubstr("Number Density is nonpositive!"))
  );
}

TEST(SourcesFactory, BulkVelocityMustBe3D) {
  const std::string model_string(
    "Constant:\n"
    "  Number Density: 1.\n"
    "  Bulk Velocity: [1., 2.]\n"
    "  Temperature: 300\n");

  const YAML::Node model_node = YAML::Load(model_string);
  const YAML::Node state_node = model_node["Constant"];

  EXPECT_THAT(
    [&]() {buildSourceStateParametersFromYAML(state_node);},
    ThrowsMessage<std::logic_error>(HasSubstr("Bulk Velocity does not have length 3!"))
  );
}

TEST(SourcesFactory, TemperatureMustBeNonNegative) {
  const std::string model_string(
    "Constant:\n"
    "  Number Density: 1.\n"
    "  Temperature: -1.\n");

  const YAML::Node model_node = YAML::Load(model_string);
  const YAML::Node state_node = model_node["Constant"];

  EXPECT_THAT(
    [&]() {buildSourceStateParametersFromYAML(state_node);},
    ThrowsMessage<std::logic_error>(HasSubstr("Temperature is negative!"))
  );
}

TEST(SourcesFactory, NumberOfParticlesMustBePositiveIfPresent) {
  const std::string model_string(
    "Initial Conditions:\n"
    "  - Species: [electron]\n"
    "    Number of Macroparticles per Species: 0\n");

  const YAML::Node model_node = YAML::Load(model_string);
  const YAML::Node initial_conditions_node = model_node["Initial Conditions"];

  EXPECT_THAT(
    [&]() {buildListOfSourceParametersFromYAML(initial_conditions_node, species_map);},
    ThrowsMessage<std::logic_error>(HasSubstr("Number of Macroparticles per Species is nonpositive!"))
  );
}

TEST(SourcesFactory, OneInitialConditionWithOneSpeciesGivesListOfOneSourceParameters) {
  constexpr double number_density = 1641356.3;
  constexpr double temperature = 300.154;
  const mfem::Vector bulk_velocity{29.3, 81.4, 95.7};
  constexpr int num_particles_per_species = 8246;
  const std::string model_string(
    "Initial Conditions:\n"
    "  - Species: [" + electron_species_name + "]\n"
    "    Number of Macroparticles per Species: " + std::to_string(num_particles_per_species) + "\n"
    "    Constant:\n"
    "      Number Density: " + std::to_string(number_density) + "\n"
    "      Temperature: " + std::to_string(temperature) + "\n"
    "      Bulk Velocity: [" + std::to_string(bulk_velocity[0]) + ", " + std::to_string(bulk_velocity[1]) + ", " + 
      std::to_string(bulk_velocity[2]) + "]\n"
  );

  const YAML::Node model_node = YAML::Load(model_string);
  const YAML::Node initial_conditions_node = model_node["Initial Conditions"];

  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters = buildListOfSourceParametersFromYAML(
    initial_conditions_node, species_map);

  EXPECT_EQ(1, list_of_parameters.size());

  const ConstantSourceParameters& parameters = dynamic_cast<const ConstantSourceParameters&>(*list_of_parameters[0]);
  EXPECT_EQ(electron_species, parameters.species);
  EXPECT_EQ(number_density, parameters.constant_state.number_density);
  EXPECT_EQ(temperature, parameters.constant_state.temperature);
  for (int i = 0; i < bulk_velocity.Size(); ++i){
    EXPECT_EQ(bulk_velocity[i], parameters.constant_state.bulk_velocity[i]);
  }
  EXPECT_EQ(num_particles_per_species, parameters.num_particles);
}

TEST(SourcesFactory, OneInitialConditionWithTwoSpeciesGivesListOfTwoSourceParameters) {
  constexpr double number_density = 2468.54;
  constexpr double temperature = 20.54;
  const mfem::Vector bulk_velocity{34.5, 72.8, 47.2};
  constexpr int num_particles_per_species = 2576;
  const std::string model_string(
    "Initial Conditions:\n"
    "  - Species: [" + electron_species_name + ", " + proton_species_name + "]\n"
    "    Number of Macroparticles per Species: " + std::to_string(num_particles_per_species) + "\n"
    "    Constant:\n"
    "      Number Density: " + std::to_string(number_density) + "\n"
    "      Temperature: " + std::to_string(temperature) + "\n"
    "      Bulk Velocity: [" + std::to_string(bulk_velocity[0]) + ", " + std::to_string(bulk_velocity[1]) + ", " + 
      std::to_string(bulk_velocity[2]) + "]\n"
  );

  const YAML::Node model_node = YAML::Load(model_string);
  const YAML::Node initial_conditions_node = model_node["Initial Conditions"];

  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters = buildListOfSourceParametersFromYAML(
    initial_conditions_node, species_map);

  EXPECT_EQ(2, list_of_parameters.size());

  for (const std::unique_ptr<SourceParameters>& parameters : list_of_parameters) {
    auto constant_parameters = dynamic_cast<const ConstantSourceParameters&>(*parameters);

    EXPECT_EQ(number_density, constant_parameters.constant_state.number_density);
    EXPECT_EQ(temperature, constant_parameters.constant_state.temperature);
    for (int i = 0; i < bulk_velocity.Size(); ++i){
      EXPECT_EQ(bulk_velocity[i], constant_parameters.constant_state.bulk_velocity[i]);
    }
    EXPECT_EQ(num_particles_per_species, parameters->num_particles);
  }

  // one species should be electron and the other should be proton
  const std::unique_ptr<SourceParameters>& parameters_0 = list_of_parameters[0];
  const std::unique_ptr<SourceParameters>& parameters_1 = list_of_parameters[1];

  EXPECT_TRUE(parameters_0->species == electron_species or parameters_0->species == proton_species);
  EXPECT_TRUE(parameters_1->species == electron_species or parameters_1->species == proton_species);
  EXPECT_NE(parameters_0->species, parameters_1->species);
}

TEST(SourcesFactory, TwoInitialConditionsWithOneSpeciesGivesListOfTwoSourceParameters) {
  constexpr double number_density_e = 2468.54;
  constexpr double temperature_e = 20.54;
  const mfem::Vector bulk_velocity_e{34.5, 72.8, 47.2};
  constexpr int num_particles_per_species_e = 2576;

  constexpr double number_density_p = 6273.54;
  constexpr double temperature_p = 432.1;
  const mfem::Vector bulk_velocity_p{82.5, 543.5, 274.5};
  constexpr int num_particles_per_species_p = 9845;

  const std::string model_string(
    "Initial Conditions:\n"
    "  - Species: [" + electron_species_name + "]\n"
    "    Number of Macroparticles per Species: " + std::to_string(num_particles_per_species_e) + "\n"
    "    Constant:\n"
    "      Number Density: " + std::to_string(number_density_e) + "\n"
    "      Temperature: " + std::to_string(temperature_e) + "\n"
    "      Bulk Velocity: [" + std::to_string(bulk_velocity_e[0]) + ", " + std::to_string(bulk_velocity_e[1]) + ", " +
      std::to_string(bulk_velocity_e[2]) + "]\n"
    "  - Species: [" + proton_species_name + "]\n"
    "    Number of Macroparticles per Species: " + std::to_string(num_particles_per_species_p) + "\n"
    "    Constant:\n"
    "      Number Density: " + std::to_string(number_density_p) + "\n"
    "      Temperature: " + std::to_string(temperature_p) + "\n"
    "      Bulk Velocity: [" + std::to_string(bulk_velocity_p[0]) + ", " + std::to_string(bulk_velocity_p[1]) + ", " + 
      std::to_string(bulk_velocity_p[2]) + "]\n"
  );

  const YAML::Node model_node = YAML::Load(model_string);
  const YAML::Node initial_conditions_node = model_node["Initial Conditions"];

  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters = buildListOfSourceParametersFromYAML(
    initial_conditions_node, species_map);

  EXPECT_EQ(2, list_of_parameters.size());

  for (const std::unique_ptr<SourceParameters>& parameters : list_of_parameters) {
    auto constant_parameters = dynamic_cast<const ConstantSourceParameters&>(*parameters);
    if (parameters->species == electron_species) {
      EXPECT_EQ(number_density_e, constant_parameters.constant_state.number_density);
      EXPECT_EQ(temperature_e, constant_parameters.constant_state.temperature);
      for (int i = 0; i < bulk_velocity_e.Size(); ++i){
        EXPECT_EQ(bulk_velocity_e[i], constant_parameters.constant_state.bulk_velocity[i]);
      }
      EXPECT_EQ(num_particles_per_species_e, parameters->num_particles);
    } else {
      EXPECT_EQ(number_density_p, constant_parameters.constant_state.number_density);
      EXPECT_EQ(temperature_p, constant_parameters.constant_state.temperature);
      for (int i = 0; i < bulk_velocity_p.Size(); ++i){
        EXPECT_EQ(bulk_velocity_p[i], constant_parameters.constant_state.bulk_velocity[i]);
      }
      EXPECT_EQ(num_particles_per_species_p, parameters->num_particles);
    }
  }

  // one species should be electron and the other should be proton
  const std::unique_ptr<SourceParameters>& parameters_0 = list_of_parameters[0];
  const std::unique_ptr<SourceParameters>& parameters_1 = list_of_parameters[1];

  EXPECT_TRUE(parameters_0->species == electron_species or parameters_0->species == proton_species);
  EXPECT_TRUE(parameters_1->species == electron_species or parameters_1->species == proton_species);
  EXPECT_NE(parameters_0->species, parameters_1->species);
}

TEST(SourcesFactory, SodInitialConditionsGivesBackCorrectParameters) {
  constexpr int num_particles_per_species = 572;

  constexpr double discontinuity_location = 0.64;

  constexpr double number_density_l = 2468.54;
  constexpr double temperature_l = 20.54;
  const mfem::Vector bulk_velocity_l{34.5, 72.8, 47.2};

  constexpr double number_density_r = 58465.54;
  constexpr double temperature_r = 320.2;
  const mfem::Vector bulk_velocity_r{-25.1, 57.4, -94.56};

  const std::string model_string(
    "Initial Conditions:\n"
    "  - Species: [" + electron_species_name + "]\n"
    "    Number of Macroparticles per Species: " + std::to_string(num_particles_per_species) + "\n"
    "    Sod:\n"
    "      Discontinuity Location: " + std::to_string(discontinuity_location) + "\n"
    "      Left State:\n"
    "        Number Density: " + std::to_string(number_density_l) + "\n"
    "        Temperature: " + std::to_string(temperature_l) + "\n"
    "        Bulk Velocity: [" + std::to_string(bulk_velocity_l[0]) + ", " + std::to_string(bulk_velocity_l[1]) + ", " +
      std::to_string(bulk_velocity_l[2]) + "]\n"
    "      Right State:\n"
    "        Number Density: " + std::to_string(number_density_r) + "\n"
    "        Temperature: " + std::to_string(temperature_r) + "\n"
    "        Bulk Velocity: [" + std::to_string(bulk_velocity_r[0]) + ", " + std::to_string(bulk_velocity_r[1]) + ", " +
      std::to_string(bulk_velocity_r[2]) + "]\n"
  );

  const YAML::Node model_node = YAML::Load(model_string);
  const YAML::Node initial_conditions_node = model_node["Initial Conditions"];

  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters = buildListOfSourceParametersFromYAML(
    initial_conditions_node, species_map);

  EXPECT_EQ(1, std::ssize(list_of_parameters));

  auto parameters = dynamic_cast<const SodSourceParameters&>(*list_of_parameters[0]);

  EXPECT_EQ(electron_species, parameters.species);
  EXPECT_EQ(num_particles_per_species, parameters.num_particles);
  EXPECT_EQ(discontinuity_location, parameters.discontinuity_location);

  EXPECT_EQ(number_density_l, parameters.left_state.number_density);
  EXPECT_EQ(temperature_l, parameters.left_state.temperature);
  for(int i = 0; i < bulk_velocity_l.Size(); ++i){
    EXPECT_EQ(bulk_velocity_l[i], parameters.left_state.bulk_velocity[i]);
  }

  EXPECT_EQ(number_density_r, parameters.right_state.number_density);
  EXPECT_EQ(temperature_r, parameters.right_state.temperature);
  for(int i = 0; i < bulk_velocity_r.Size(); ++i){
    EXPECT_EQ(bulk_velocity_r[i], parameters.right_state.bulk_velocity[i]);
  }
}

TEST(SourcesFactory, GaussianInitialConditionsGivesBackCorrectParameters) {
  constexpr int num_particles_per_species = 513;
  constexpr double center = 0.45;
  constexpr double standard_deviation = 0.014;

  constexpr double number_density_offset = 1e19;
  constexpr double number_density_height = 5e18;
  constexpr double pressure_offset = 0.04;
  constexpr double pressure_height = 0.005;
  const mfem::Vector bulk_velocity_offset{10.3, 15.2, 19.5};
  const mfem::Vector bulk_velocity_height{2.4, 3.1, 1.9};

  const std::string model_string(
    "Initial Conditions:\n"
    "  - Species: [" + proton_species_name + "]\n"
    "    Number of Macroparticles per Species: " + std::to_string(num_particles_per_species) + "\n"
    "    Gaussian:\n"
    "      Center: [" + std::to_string(center) + "]\n"
    "      Standard Deviation: " + std::to_string(standard_deviation) + "\n"
    "      Offsets:\n"
    "        Number Density: " + std::to_string(number_density_offset) + "\n"
    "        Pressure: " + std::to_string(pressure_offset) + "\n"
    "        Bulk Velocity: [" + std::to_string(bulk_velocity_offset[0]) + ", " + std::to_string(bulk_velocity_offset[1]) + ", " +
      std::to_string(bulk_velocity_offset[2]) + "]\n"
    "      Heights:\n"
    "        Number Density: " + std::to_string(number_density_height) + "\n"
    "        Pressure: " + std::to_string(pressure_height) + "\n"
    "        Bulk Velocity: [" + std::to_string(bulk_velocity_height[0]) + ", " + std::to_string(bulk_velocity_height[1]) + ", " +
      std::to_string(bulk_velocity_height[2]) + "]\n"
  );

  const YAML::Node model_node = YAML::Load(model_string);
  const YAML::Node initial_conditions_node = model_node["Initial Conditions"];

  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters = buildListOfSourceParametersFromYAML(
    initial_conditions_node, species_map);

  EXPECT_EQ(1, std::ssize(list_of_parameters));

  auto parameters = dynamic_cast<const GaussianSourceParameters&>(*list_of_parameters[0]);

  EXPECT_EQ(proton_species, parameters.species);
  EXPECT_EQ(num_particles_per_species, parameters.num_particles);
  EXPECT_EQ(1, parameters.center.Size());
  EXPECT_EQ(center, parameters.center[0]);

  EXPECT_EQ(number_density_offset, parameters.offsets.number_density);
  EXPECT_EQ(0, parameters.offsets.temperature);
  for (int i = 0; i < bulk_velocity_offset.Size(); ++i) {
    EXPECT_EQ(bulk_velocity_offset[i], parameters.offsets.bulk_velocity[i]);
  }
  EXPECT_EQ(pressure_offset, parameters.pressure_offset);

  EXPECT_EQ(number_density_height, parameters.heights.number_density);
  EXPECT_EQ(0, parameters.heights.temperature);
  for (int i = 0; i < bulk_velocity_height.Size(); ++i) {
    EXPECT_EQ(bulk_velocity_height[i], parameters.heights.bulk_velocity[i]);
  }
  EXPECT_EQ(pressure_height, parameters.pressure_height);
}

}