#include <H5Fpublic.h>
#include <H5Ppublic.h>
#include <H5Spublic.h>
#include <H5Tpublic.h>
#include <libmfpic/DumpParticles.hpp>
#include <libmfpic/ParticleContainer.hpp>

#include <hdf5.h>

namespace mfpic {

static bool file_is_created = false;

void dumpParticles(const ParticleContainer& particles, double simulation_time, const std::string filename) {
  hid_t file;
  if (file_is_created) {
    file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  } else {
    file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    file_is_created = true;
  }

  std::vector<double> x, y, z, vx, vy, vz, weight;
  std::vector<int> element;
  for (const Particle& particle : particles) {
    if (particle.is_alive) {
      x.push_back(particle.position[0]);
      y.push_back(particle.position[1]);
      z.push_back(particle.position[2]);
      vx.push_back(particle.velocity[0]);
      vy.push_back(particle.velocity[1]);
      vz.push_back(particle.velocity[2]);
      weight.push_back(particle.weight);
      element.push_back(particle.element);
    }
  }

  H5G_info_t top_level_group_info;
  H5Gget_info_by_name(file, "/", &top_level_group_info, H5P_DEFAULT);
  const int step_to_write = top_level_group_info.nlinks;

  std::string step_name = "Step#" + std::to_string(step_to_write);
  hid_t step_group = H5Gcreate(file, step_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  constexpr hsize_t one = 1;
  hid_t timevalue_dataspace = H5Screate_simple(1, &one, NULL);
  hid_t timevalue_attribute = H5Acreate(step_group, "TimeValue", H5T_NATIVE_DOUBLE, timevalue_dataspace, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(timevalue_attribute, H5T_NATIVE_DOUBLE, &simulation_time);
  H5Aclose(timevalue_attribute);
  H5Sclose(timevalue_dataspace);

  const hsize_t num_particles = x.size();
  hid_t dataspace = H5Screate_simple(1, &num_particles, NULL);
  hid_t dataset = H5Dcreate(step_group, "x", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, x.data());
  H5Dclose(dataset);
  dataset = H5Dcreate(step_group, "y", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, y.data());
  H5Dclose(dataset);
  dataset = H5Dcreate(step_group, "z", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, z.data());
  H5Dclose(dataset);
  dataset = H5Dcreate(step_group, "vx", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vx.data());
  H5Dclose(dataset);
  dataset = H5Dcreate(step_group, "vy", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vy.data());
  H5Dclose(dataset);
  dataset = H5Dcreate(step_group, "vz", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vz.data());
  H5Dclose(dataset);
  dataset = H5Dcreate(step_group, "weight", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, weight.data());
  H5Dclose(dataset);
  dataset = H5Dcreate(step_group, "element", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, element.data());
  H5Dclose(dataset);

  H5Sclose(dataspace);
  H5Gclose(step_group);
  H5Fclose(file);
}

} // namespace mfpic
