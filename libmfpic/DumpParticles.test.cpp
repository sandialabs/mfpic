#include <H5Tpublic.h>
#include <libmfpic/DumpParticles.hpp>
#include <libmfpic/ParticleContainer.hpp>

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
    });
  }

  return std::make_pair(timevalue, particles);
}

TEST(DumpParticles, ReadParticlesMatcDumpedParticles) {
  const mfem::Vector position({1.0, 2.0, 3.0});
  const mfem::Vector velocity({4.0, 5.0, 6.0});
  constexpr int element = 505;
  constexpr double weight = 1.0e9;
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
      EXPECT_EQ(element, particle.element);
      for (int dim = 0; dim < 3; dim++) {
        EXPECT_DOUBLE_EQ(position[dim], particle.position[dim]);
        EXPECT_DOUBLE_EQ(velocity[dim], particle.velocity[dim]);
      }
    }
  }
}

} // namespace
