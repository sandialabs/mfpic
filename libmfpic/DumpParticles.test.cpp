#include <H5Tpublic.h>
#include <libmfpic/DGEulerAssembly.hpp>
#include <libmfpic/DGEulerOperations.hpp>
#include <libmfpic/DGEulerInitialConditionsFactory.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/DumpParticles.hpp>
#include <libmfpic/LoadParticles.hpp>
#include <libmfpic/ParticleContainer.hpp>
#include <libmfpic/ParticleOperations.hpp>
#include <libmfpic/ReflectingParticleBoundary.hpp>
#include <libmfpic/SourcesFactory.hpp>
#include <libmfpic/Species.hpp>

#include <gtest/gtest.h>
#include <hdf5.h>

namespace {

using namespace mfpic;

constexpr char filename[] = "particles.h5part";

std::pair<double, ParticleContainer> readTimeValueAndParticlesFromStep(int step) {
  hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  const std::string step_name = "Step#" + std::to_string(step);
  hid_t step_group = H5Gopen(file, step_name.c_str(), H5P_DEFAULT);

  hid_t timevalue_attribute = H5Aopen(step_group, "TimeValue", H5P_DEFAULT);
  double timevalue;
  H5Aread(timevalue_attribute, H5T_NATIVE_DOUBLE, &timevalue);
  H5Aclose(timevalue_attribute);

  hid_t dataset = H5Dopen(step_group, "x", H5P_DEFAULT);
  hid_t dataspace = H5Dget_space(dataset);
  const int num_particles = H5Sget_simple_extent_npoints(dataspace);
  H5Sclose(dataspace);

  std::vector<double> x(num_particles);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, x.data());
  H5Dclose(dataset);
  dataset = H5Dopen(step_group, "y", H5P_DEFAULT);
  std::vector<double> y(num_particles);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, y.data());
  H5Dclose(dataset);
  dataset = H5Dopen(step_group, "z", H5P_DEFAULT);
  std::vector<double> z(num_particles);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, z.data());
  H5Dclose(dataset);
  dataset = H5Dopen(step_group, "vx", H5P_DEFAULT);
  std::vector<double> vx(num_particles);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vx.data());
  H5Dclose(dataset);
  dataset = H5Dopen(step_group, "vy", H5P_DEFAULT);
  std::vector<double> vy(num_particles);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vy.data());
  H5Dclose(dataset);
  dataset = H5Dopen(step_group, "vz", H5P_DEFAULT);
  std::vector<double> vz(num_particles);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vz.data());
  H5Dclose(dataset);
  dataset = H5Dopen(step_group, "weight", H5P_DEFAULT);
  std::vector<double> weight(num_particles);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, weight.data());
  H5Dclose(dataset);
  dataset = H5Dopen(step_group, "particle_distribution_function_value", H5P_DEFAULT);
  std::vector<double> particle_distribution_function_value(num_particles);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, particle_distribution_function_value.data());
  H5Dclose(dataset);
  dataset = H5Dopen(step_group, "element", H5P_DEFAULT);
  std::vector<int> element(num_particles);
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, element.data());
  H5Dclose(dataset);

  H5Gclose(step_group);
  H5Fclose(file);

  ParticleContainer particles;
  for (int i = 0; i < num_particles; i++) {
    particles.addParticle(Particle{
      .position = mfem::Vector({x[i], y[i], z[i]}),
      .velocity = mfem::Vector({vx[i], vy[i], vz[i]}),
      .element = element[i],
      .weight = weight[i],
      .is_alive = true,
      .particle_distribution_function_value = particle_distribution_function_value[i],
    });
  }

  return std::make_pair(timevalue, particles);
}

TEST(DumpParticles, ReadParticlesMatchDumpedParticles) {
  const mfem::Vector position({1.0, 2.0, 3.0});
  const mfem::Vector velocity({4.0, 5.0, 6.0});
  constexpr int element = 505;
  constexpr double weight = 1.0e9;
  constexpr double particle_distribution_function_value = 21492e9;
  ParticleContainer particles_to_dump;
  constexpr int num_timesteps = 3;
  constexpr double dt = 1.0e-12;
  for (int i = 0; i < num_timesteps; i++) {
    const double simulation_time = dt * i;
    dumpParticles(particles_to_dump, simulation_time);
    particles_to_dump.addParticle(Particle{
      .position = position,
      .velocity = velocity,
      .element = element,
      .weight = weight,
      .particle_distribution_function_value = particle_distribution_function_value
    });
  }

  for (int i = 0; i < num_timesteps; i++) {
    auto [simulation_time, particles_from_step] = readTimeValueAndParticlesFromStep(i);
    const int expected_num_particles = i;
    ASSERT_EQ(expected_num_particles, particles_from_step.numParticles());
    const double expected_time = i * dt;
    EXPECT_DOUBLE_EQ(expected_time, simulation_time);
    for (const Particle& particle : particles_from_step) {
      EXPECT_DOUBLE_EQ(weight, particle.weight);
      EXPECT_DOUBLE_EQ(particle_distribution_function_value, particle.particle_distribution_function_value);
      EXPECT_EQ(element, particle.element);
      for (int dim = 0; dim < 3; dim++) {
        EXPECT_DOUBLE_EQ(position[dim], particle.position[dim]);
        EXPECT_DOUBLE_EQ(velocity[dim], particle.velocity[dim]);
      }
    }
  }
}


static void removeIfExists(const std::string& path)
{
  std::remove(path.c_str());
}

static std::vector<std::string> readLines(const std::string& path)
{
  std::ifstream in(path);
  EXPECT_TRUE((bool)in) << "Failed to open: " << path;

  std::vector<std::string> lines;
  std::string line;
  while (std::getline(in, line)) lines.push_back(line);
  return lines;
}

static std::vector<std::string> splitCsv(const std::string& line)
{
  std::vector<std::string> toks;
  std::stringstream ss(line);
  std::string tok;
  while (std::getline(ss, tok, ',')) toks.push_back(tok);
  return toks;
}

static int stringToInt(const std::string& s) { return std::stoi(s); }
static double stringToDouble(const std::string& s) { return std::stod(s); }

constexpr Species default_species{.charge = 1.0, .mass = 1.0, .id = 0};
const std::vector<std::shared_ptr<ParticleBoundaryFactory>> empty_particle_boundary_factory_list;
const std::shared_ptr<ParticleBoundaryFactory> default_reflecting_particle_boundary_factory
  = std::make_shared<ReflectingParticleBoundaryFactory>();
constexpr int one_species = 1, two_species = 2;

TEST(DumpParticles, CSVCorrectForOneSpecies)
{
  const int num_elems = 5;
  std::shared_ptr<mfem::Mesh> mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elems, .234));
  constexpr int order = 1;
  Discretization discretization(mesh.get(),order);

  const mfem::Vector nominal_bulk_velocity({300.0,400.0,500.0});
  constexpr double temperature = 11600.0;
  constexpr double number_density = 1.0e18;
  const SourceStateParameters source_state_parameters{
    .number_density = number_density,
    .bulk_velocity = nominal_bulk_velocity,
    .temperature = temperature,
  };
  constexpr int num_particles = 100;
  constexpr int num_species = 1; 

  std::mt19937 generator;

  ParticleContainer particles = loadParticles(
    ConstantSourceParameters(default_species, source_state_parameters, num_particles),
    generator,
    mesh  
  );

  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );

  const std::string prefix = "csv_test_one";

  for (int s = 0; s < num_species; ++s)
    removeIfExists(prefix + "_species_" + std::to_string(s) + ".csv");

  const int steps = 10;
  const double dt = 0.25;

  mfem::DenseMatrix& computed_number_density = particle_operations.getNumberDensity(particles);
  mfem::DenseTensor& computed_bulk_velocity  = particle_operations.getBulkVelocity(particles);
  mfem::DenseMatrix& computed_temperature  = particle_operations.getTemperature(particles);

  for (int step = 0; step < steps; ++step)
    dumpParticleMoments(particle_operations,particles, prefix, step, dt * step);

  for (int s = 0; s < num_species; ++s) {
    const auto file = prefix + "_species_" + std::to_string(s) + ".csv";
    const auto lines = readLines(file);
    ASSERT_EQ(lines.size(), static_cast<size_t>(1 + num_elems*steps));

    for (int step = 0; step < steps; ++step) {
      for (int e = 0; e < num_elems; ++e) {
        const std::string& line = lines[1 + step*num_elems + e];
        auto toks = splitCsv(line);
        ASSERT_EQ(toks.size(), static_cast<size_t>(11));

        EXPECT_EQ(stringToInt(toks[0]), step);
        EXPECT_DOUBLE_EQ(stringToDouble(toks[1]), dt*step);
        EXPECT_EQ(stringToInt(toks[2]), e);

        const double csv_number_density = stringToDouble(toks[6]);
        const double csv_temperature = stringToDouble(toks[7]);
        const double csv_bulk_velocity_0 = stringToDouble(toks[8]);
        const double csv_bulk_velocity_1 = stringToDouble(toks[9]);
        const double csv_bulk_velocity_2 = stringToDouble(toks[10]);

        EXPECT_NEAR(csv_number_density, computed_number_density(e, s), 1e-12);
        EXPECT_NEAR(csv_temperature,  computed_temperature(e, s),  1e-12);

        const mfem::Vector computed_bulk_velocity_in_element(computed_bulk_velocity(s).GetColumn(e), 3);
        EXPECT_NEAR(csv_bulk_velocity_0, computed_bulk_velocity_in_element(0), 1e-12);
        EXPECT_NEAR(csv_bulk_velocity_1, computed_bulk_velocity_in_element(1), 1e-12);
        EXPECT_NEAR(csv_bulk_velocity_2, computed_bulk_velocity_in_element(2), 1e-12);
      }
    }
  }
};

TEST(DumpParticles, CSVCorrectForTwoSpecies)
{

  constexpr Species species_1{.charge = 1.0, .mass = 1.0, .id = 0};
  constexpr Species species_2{.charge = 2.0, .mass = 2.0, .id = 1};
  const int num_elems = 5;
  std::shared_ptr<mfem::Mesh> mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elems, .234));
  constexpr int order = 1;
  Discretization discretization(mesh.get(),order);

  const mfem::Vector nominal_bulk_velocity({300.0,400.0,500.0});
  constexpr double temperature = 11600.0;
  constexpr double number_density = 1.0e18;
  const SourceStateParameters source_state_parameters{
    .number_density = number_density,
    .bulk_velocity = nominal_bulk_velocity,
    .temperature = temperature,
  };
  constexpr int num_particles = 100;
  constexpr int num_species = 2; 

  std::mt19937 generator;

  ParticleContainer particles = loadParticles(
    ConstantSourceParameters(species_1, source_state_parameters, num_particles),
    generator,
    mesh  
  );

  particles.addParticles(loadParticles(
    ConstantSourceParameters(species_2, source_state_parameters, num_particles),
    generator,
    mesh
  ));

  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    two_species
  );

  const std::string prefix = "csv_test_two";

  for (int s = 0; s < num_species; ++s)
    removeIfExists(prefix + "_species_" + std::to_string(s) + ".csv");

  const int steps = 10;
  const double dt = 0.25;

  mfem::DenseMatrix& computed_number_density = particle_operations.getNumberDensity(particles);
  mfem::DenseTensor& computed_bulk_velocity  = particle_operations.getBulkVelocity(particles);
  mfem::DenseMatrix& computed_temperature  = particle_operations.getTemperature(particles);

  for (int step = 0; step < steps; ++step)
    dumpParticleMoments(particle_operations,particles, prefix, step, dt * step);

  for (int s = 0; s < num_species; ++s) {
    const auto file = prefix + "_species_" + std::to_string(s) + ".csv";
    const auto lines = readLines(file);
    ASSERT_EQ(lines.size(), static_cast<size_t>(1 + num_elems*steps));

    for (int step = 0; step < steps; ++step) {
      for (int e = 0; e < num_elems; ++e) {
        const std::string& line = lines[1 + step*num_elems + e];
        auto toks = splitCsv(line);
        ASSERT_EQ(toks.size(), static_cast<size_t>(11));

        EXPECT_EQ(stringToInt(toks[0]), step);
        EXPECT_DOUBLE_EQ(stringToDouble(toks[1]), dt*step);
        EXPECT_EQ(stringToInt(toks[2]), e);

        const double csv_number_density = stringToDouble(toks[6]);
        const double csv_temperature = stringToDouble(toks[7]);
        const double csv_bulk_velocity_0 = stringToDouble(toks[8]);
        const double csv_bulk_velocity_1 = stringToDouble(toks[9]);
        const double csv_bulk_velocity_2 = stringToDouble(toks[10]);

        EXPECT_NEAR(csv_number_density, computed_number_density(e, s), 1e-12);
        EXPECT_NEAR(csv_temperature,  computed_temperature(e, s),  1e-12);

        const mfem::Vector computed_bulk_velocity_in_element(computed_bulk_velocity(s).GetColumn(e), 3);
        EXPECT_NEAR(csv_bulk_velocity_0, computed_bulk_velocity_in_element(0), 1e-12);
        EXPECT_NEAR(csv_bulk_velocity_1, computed_bulk_velocity_in_element(1), 1e-12);
        EXPECT_NEAR(csv_bulk_velocity_2, computed_bulk_velocity_in_element(2), 1e-12);
      }
    }
  }
};

} // namespace
